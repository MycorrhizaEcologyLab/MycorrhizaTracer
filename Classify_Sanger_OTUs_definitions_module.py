################################################
#imports
import os
import glob
import sys
import pprint
import re
from multiprocessing import Pool
import subprocess
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord	
from statistics import mode
from statistics import median
import requests
from Bio.Align import PairwiseAligner
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

################################################################################################
#Definitions:

#####SETUP#######################################################

def make_output_dirs(outdir):
	print("\nBegin by building an output file heirarchy if none already exists.", file=sys.stderr)
	#needs to make the appropriate output directory structure
	
	dirs_to_make = [
		outdir, 
		outdir + "/01_consensuses/", 
		outdir + "/02_classify_seqs/", 
		outdir + "/03_cluster_seqs/",
		outdir + "/04_FUNGuild/"
	]
	
	for dir_to_make in dirs_to_make:
		if not os.path.exists(dir_to_make):
			try:
				os.makedirs(dir_to_make)
				print(f"Created directory: {dir_to_make}", file=sys.stderr)
			except OSError as e:
				print(f"Error creating directory {dir_to_make}: {e}", file=sys.stderr)
				sys.exit(1) 
	print("\tdone.", file=sys.stderr)

def add_in_dirDepth(SampleID, gene, args):
	"""Returns the subdirectory path (after 01_consensuses/ and before the filename) based on dirDepth logic."""
	"""this helps build path names in many places and keep track of intermediary files"""
	fileName = f"{SampleID}_{gene}.consensus.fq"
	if hasattr(args, 'dirDepth') and args.dirDepth:
		bits = fileName.split("_")
		dirDepth = min(len(bits)-1, args.dirDepth)
		subdir = os.path.join(*bits[:dirDepth]) if dirDepth > 0 else ""
		return subdir
	return ""

def load_metadata(args):
	print("Reading in Metadata file: ", args.metadata, sep="", file=sys.stderr)
	resequenced_counter = linecounter=0
	
	gene_list = ["ITS", "ITS2", "RBCL"]
	if args.onlyITS:
		gene_list = ["ITS"]
	
	linecounter = 0
	MetaDict = {}

	with open (args.metadata, "r", encoding='utf-8-sig') as FH:
		for line in FH:
			#skip the header line
			#input(line)
			if line.strip().startswith("Sample_ID"):
				#input("HEADER")
				continue
			linecounter += 1
			if linecounter % 10 ==0:
				print("	Samples processed: ", linecounter, end = "\r", file=sys.stderr, sep="")
			if linecounter > args.firstNsamples and args.firstNsamples > 0:
				print(f"Processed {args.firstNsamples} samples. Stopping early as requested.", file=sys.stderr)
				input("Press Enter to continue...")  # Debugging pause
				break
			line=line.strip("\n")
			bits = line.split(",")
			#skip empty lines (typically ones with just commas)
			if len(bits) < 9: #implies that the line is empty - only a bunch of commas. could also be solved by removing empty lines from csv file
				continue
	
			#Sample_ID, SalvageGroup, directory, ITS_F, ITS_R, ITS2_F, ITS2_R, RBCL_F, RBCL_R
			
			Sample_ID =		bits[0]
			SalvageGroup =	bits[1]
			directory =		bits[2]
			ITS_F =			bits[3]
			ITS_R =			bits[4]
			ITS2_F = 		bits[5]
			ITS2_R = 		bits[6]
			RBCL_F = 		bits[7]
			RBCL_R = 		bits[8]
			
			#add in things to metadict: this is where it's all stored - data, metadata, completed steps, everything. 
			if Sample_ID not in MetaDict:
				MetaDict[Sample_ID] = {
					"Metadata":{
						"SalvageGroup":SalvageGroup,
						"raw_seqs_dir":directory,
					},	
					"ITS":{
						"Seqs":{
							"ab1s":{
								"records":[],
								"ab1_files":[],
								"fq_files":[],
								"seqLens":[],
								"avQuals":[],
							},
							"Cons":{
								"record":None,
								"fasta":None,
								"fastq":None,
								"seqLen":None,
								"avQual":None,
								"alignmentQual":None,
								"nSeqs":None,
							}
						},
						"Classifications":{#using a blast to the UNITE/UNITE_all/RBCL databases and then clustered
							"classification_type":None, 
							"blastfile":None,
							"taxonomy":None, 
							"best_hit":None,
							"SH_number":None
						},
						"Clustering":{
							"ITSx_trim_queue":False, #this is set to True if the sequence was passed to ITSX
							"ITSx_trimmed":False, #this is set to True if the sequences was successfully trimmed by ITSX - i.e.: if ITSX found ITS in it. 
							"Send_to_Clustering":False,
							"Clustered":False, #this is set to True if the sequence was successfully clustered OR salvaged into a cluster: so the number of trues of this + the number of trues of clustered_by_salvage should match ITSx_trimmed....
							"clustered_by_salvage":False, #this is set to True if the sequence was successfully salvaged. The total of these that are true should be less than (or equal to the number of trues of ITSx_queue - the number of trues of ITSx_trimmed)
							"cluster_data":{
								"qseqid":None, 
								"clusterID":None, 
								"centroid":None, 
								"cigar":None, 
								"pident":None, 
								"length":None, 
								"clustSize":None, 
								"strand":None
							}
						},
						"FUNGuild":{
							"Taxon":"NA",
							"Taxon_Level":"NA",
							"Trophic_Mode":"NA",
							"Guild":"NA",
							"Growth_Morphology":"NA",
							"Trait":"NA",
							"Confidence_Ranking":"NA",
							"Notes":"NA",
							"Citation_Source":"NA"
						},
						"Salvaging":{
							"attempted":False, #this is set to True if the sequence was salvaged from a low quality read
							"Successfull":False, #this is set to True if the sequence was salvaged from a low quality read: should match "clustered_by_salvage"
						}
					},
					"ITS2":{
						"Seqs":{
							"ab1s":{
								"records":[],
								"ab1_files":[],
								"fq_files":[],
								"seqLens":[],
								"avQuals":[],
							},
							"Cons":{
								"record":None,
								"fasta":None,
								"fastq":None,
								"seqLen":None,
								"avQual":None,
								"alignmentQual":None,
								"nSeqs":None,
							}
						},
						"Classifications":{#using a blast to the UNITE/UNITE_all/RBCL databases and then clustered
							"classification_type":None, #this is set to True if the sequence was salvaged from a low quality read
							"blastfile":None,
							"taxonomy":None, 
							"best_hit":None,
							"SH_number":None
						},
						"Clustering":{
							"ITSx_trim_queue":False, #this is set to True if the sequence was passed to ITSX
							"ITSx_trimmed":False, #this is set to True if the sequences was successfully trimmed by ITSX - i.e.: if ITSX found ITS in it. 
							"Send_to_Clustering":False,
							"Clustered":False, #this is set to True if the sequence was successfully clustered OR salvaged into a cluster: so the number of trues of this + the number of trues of clustered_by_salvage should match ITSx_trimmed....
							"clustered_by_salvage":False, #this is set to True if the sequence was successfully salvaged. The total of these that are true should be less than (or equal to the number of trues of ITSx_queue - the number of trues of ITSx_trimmed)
							"cluster_data":{
								"qseqid":None, 
								"clusterID":None, 
								"centroid":None, 
								"cigar":None, 
								"pident":None, 
								"length":None, 
								"clustSize":None, 
								"strand":None
							}
						}
					},
					"RBCL":{
						"Seqs":{
							"ab1s":{
								"records":[],
								"ab1_files":[],
								"fq_files":[],
								"seqLens":[],
								"avQuals":[],
							},
							"Cons":{
								"record":None,
								"fasta":None,
								"fastq":None,
								"seqLen":None,
								"avQual":None,
								"alignmentQual":None,
								"nSeqs":None,
							}
						},
						"Classifications":{#using a blast to the UNITE/UNITE_all/RBCL databases and then clustered
							"classification_type":None, #this is set to True if the sequence was salvaged from a low quality read
							"blastfile":None,
							"taxonomy":None, 
							"best_hit":None,
							"SH_number":None
						},
						"Clustering":{
							#"ITSx_trim_queue":False, #this is set to True if the sequence was passed to ITSX
							#"ITSx_trimmed":False, #this is set to True if the sequences was successfully trimmed by ITSX - i.e.: if ITSX found ITS in it. 
							"Send_to_Clustering":False,
							"Clustered":False, #this is set to True if the sequence was successfully clustered OR salvaged into a cluster: so the number of trues of this + the number of trues of clustered_by_salvage should match ITSx_trimmed....
							"clustered_by_salvage":False, #this is set to True if the sequence was successfully salvaged. The total of these that are true should be less than (or equal to the number of trues of ITSx_queue - the number of trues of ITSx_trimmed)
							"cluster_data":{
								"qseqid":None, 
								"clusterID":None, 
								"centroid":None, 
								"cigar":None, 
								"pident":None, 
								"length":None, 
								"clustSize":None, 
								"strand":None
							}
						}
					}
				}
			flag = False
			
			for gene in gene_list:
				#check that at least one exists (by skipping the sample they are both absent.) and then add the ab1 filepaths to the metadata dictionary:
				if not os.path.exists(os.path.join(args.prefix, directory, eval(f"{gene}_F"))) and not os.path.exists(os.path.join(args.prefix, directory, eval(f"{gene}_R"))):
					print(os.path.join(args.prefix, directory, eval(f"{gene}_F")), "does not exist", file=sys.stderr)
					flag = True
					continue
				#elif not os.path.exists(os.path.join(args.prefix, directory, eval(f"{gene}_R"))):
				#	pprint(os.path.join(args.prefix, directory, eval(f"{gene}_R")), "does not exist", file=sys.stderr)
				#	flag = True
				#	continue
				else:
					#add the ab1 filepaths to the metadata dictionary: #this might break if there is only a forward or reverse read for a sample, need to make it more resilient to that.
					#print(os.path.join(args.prefix, directory, eval(f"{gene}_F")))
					forward_read_present = os.path.isfile(os.path.join(args.prefix, directory, eval(f"{gene}_F")))
					reverse_read_present = os.path.isfile(os.path.join(args.prefix, directory, eval(f"{gene}_R")))
					#print("forward read present?", os.path.join(args.prefix, directory, eval(f"{gene}_F")))
					#input(forward_read_present)
					#print("reverse read present?", os.path.join(args.prefix, directory, eval(f"{gene}_R")))
					#input(reverse_read_present)
					if forward_read_present:
						MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["ab1_files"].append(os.path.join(args.prefix, directory, eval(f"{gene}_F")))
					if reverse_read_present:
						MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["ab1_files"].append(os.path.join(args.prefix, directory, eval(f"{gene}_R")))
				
				#add the fastq filepaths to the metadata dictionary:
				MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"].append(os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_"+gene+".Forward.fq"))
				MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"].append(os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_"+gene+".Reverse.fq"))

				#add the consensus fastq filepaths to the metadata dictionary:
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"] = os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_"+gene+".consensus.fa")
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fastq"] = os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_"+gene+".consensus.fq")
			if flag:
				print(f"Skipping Sample ID: {Sample_ID} due to missing files.", file=sys.stderr)
				continue
	print("	Samples processed: ", linecounter, end = "\r", file=sys.stderr, sep="")
	return(MetaDict)

#####CONSENSUS#######################################################
def create_fastq_from_ab1(MetaDict, args, SS_FH):
	print("\nCreating fastq files from ab1 files.", file=sys.stderr)
	
	gene_list = ["ITS", "ITS2", "RBCL"]
	if args.onlyITS:
		gene_list = ["ITS"]


	counter = 0
	for Sample_ID in MetaDict:
		counter +=1
		if counter % 10 == 0:
			print("	Samples processed: ", counter, end = "\r", file=sys.stderr, sep="")
		if args.verbose:
			print(f"Processing Sample_ID: {Sample_ID}", file=sys.stderr)
		for gene in gene_list:
			#convert ab1 files to fastq files using Bio.SeqIO
			#print(Sample_ID, gene)
			#pprint.pprint(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"])
			#input("Press Enter to continue...")
			for ab1_file,read in zip(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["ab1_files"], range(0,len(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["ab1_files"]))): 
				if not os.path.exists(os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args))):
						os.makedirs(os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args)))

				fq_filename = MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"][read]
				if args.resume and os.path.exists(fq_filename):
					if args.verbose:
						print(f"{fq_filename} already exists. Loading this file from memory.", file=sys.stderr)
					#read it in, and. update metadict, then continue to the next file
					seq_record = SeqIO.read(fq_filename, "fastq")
					MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["records"].append(seq_record)
					seqLen = len(seq_record)
					seqQual = round(sum(seq_record.letter_annotations["phred_quality"]) / len(seq_record.letter_annotations["phred_quality"]),1)
					MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["seqLens"].append(seqLen)
					MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["avQuals"].append(seqQual)
					if args.sequence_stats:
						#write the sequence stats to a file
						print(Sample_ID, gene, read, seqLen, seqQual, sep="	", file=SS_FH)
					continue
				
				if args.verbose:
					print(f"{fq_filename} needs creating. Importing and trimming this file.", file=sys.stderr)

				if read == 0:
					FR = "Forward"
				elif read == 1:
					FR = "Reverse"
				if args.troubleshooting:
					#I would like to cp ab1_file to the directory where the fastq file will be written, so that I can inspect it.			
					cp_cmd = ["cp", ab1_file, os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_"+gene+"."+FR+".ab1")]
					subprocess.run(cp_cmd, check=True, capture_output=True)
				#read the ab1 file and write it to a fastq file
				#pprint.pprint(ab1_file)
				#pprint.pprint(Sample_ID)
				#input("Press Enter to continue...")  # Debugging pause
				with open(ab1_file, "rb") as ab1_handle:
					#try to read the ab1 file in ABI format
					try:
						seq_record = SeqIO.read(ab1_handle, "abi") #this one just reads in the ab1 file without trimming by peak height
						#print(seq_record.annotations["abif_raw"]["DATA9"], file=sys.stderr) #debugging line to see the seq_record
						#input("Press Enter to continue...")  # Debugging pause
					except Exception as e:
						try:
							ab1_handle.seek(0) # Reset file pointer to the beginning
							seq_record = SeqIO.read(ab1_handle, "scf")
							print(f"Warning: Failed to read {ab1_file} as ABI, read as SCF instead. Error: {e}", file=sys.stderr)
						except Exception as e2:
							print(f"Error: Failed to read {ab1_file} as ABI or SCF. ABI error: {e}, SCF error: {e2}", file=sys.stderr)
							raise
					
					#need to qual trim the sequence
					seq_record = trim_seq_record(seq_record, args)
					seq_record = trim_seq_by_peak_height(seq_record, window_size=args.window_size, stddev_cutoff=args.stddev_cutoff) #this reads in the ab1 and trims by peak height. 
					#input("Write a primer-aligner trim function that looks for the RC of R pimer in F read, and the RC of F primer in the R Read and trims off any trailing sequence.")
					#exit()
					seq_record.id = Sample_ID
					seq_record.description = gene+" "+FR #add the gene and read number to the description
					if read % 2 == 1: #these should genearlly be reverse reads, so reverse complement them; might fail if there are more than 1 sequencing attempt and different F and R reads failed in the different attempts, but that seems super unlikely.
						rcseq_record = seq_record.reverse_complement() #now need to copy metadata from the original record to the new one
						rcseq_record.id = seq_record.id
						rcseq_record.description = seq_record.description + " revcom"
						rcseq_record.name = seq_record.name + " revcom" #not used for seqio.write
						seq_record = rcseq_record
					#add the sequence record to the metadata dictionary
					MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["records"].append(seq_record)
					seqLen = len(seq_record)
					seqQual = round(sum(seq_record.letter_annotations["phred_quality"]) / len(seq_record.letter_annotations["phred_quality"]),1)
					#add the sequence length and
					MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["seqLens"].append(seqLen)
					MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["avQuals"].append(seqQual)
					if args.sequence_stats:
						#write the sequence stats to a file
						print(Sample_ID, gene, read, seqLen, seqQual, sep="	", file=SS_FH)

					#pprint.pprint(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"])
					if not os.path.exists(os.path.dirname(fq_filename)):
						os.makedirs(os.path.dirname(fq_filename))
					with open(fq_filename, "w") as fq_handle:
						SeqIO.write(seq_record, fq_handle, "fastq")
				#if gene== "ITS":
				#	pprint.pprint(seq_record)
				#	input("Press Enter to continue...")
	#pprint.pprint(MetaDict["N692_4"]["ITS"])
	#input("1! does it have a consensus fasta/q? no quals")


	print("	Samples processed: ", counter, end = "\r", file=sys.stderr, sep="")
	return MetaDict

def make_consensus_fastq(MetaDict, args):
	print("\nMaking consensus fasta and fastq files from ab1 files.", file=sys.stderr)
	gene_list = ["ITS", "ITS2", "RBCL"]
	if args.onlyITS:
		gene_list = ["ITS"]
	counter = 0
	if args.manual_consensus_file:
		#need to read in the manual consensus file and add those sequences to the metadict instead of making new ones.
		#so read it into a dict, so that each time we try a new sample*gene we can check if it's in the manual consensus dict first.
		manual_consensus_dict = read_in_manual_consensus_file(args)
		
		# Check which manual consensus samples are NOT in the main metadata
		missing_from_metadata = [sid for sid in manual_consensus_dict if sid not in MetaDict]
		if missing_from_metadata:
			print(f"\nWarning: {len(missing_from_metadata)} Sample_IDs from manual consensus file are not found in the main metadata file:", file=sys.stderr)
			for sid in missing_from_metadata[:10]:  # Show first 10
				print(f"  - {sid}", file=sys.stderr)
			if len(missing_from_metadata) > 10:
				print(f"  ... and {len(missing_from_metadata) - 10} more", file=sys.stderr)
			print("  These may be typos or samples that were filtered out. Check for exact Sample_ID matches.", file=sys.stderr)

	for Sample_ID in MetaDict:
		counter+=1
		if counter % 10 == 0:
			print("	Samples processed: ", counter, end = "\r", file=sys.stderr, sep="")
		for gene in gene_list:
			consensus_fastq = MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fastq"]
			if args.manual_consensus_file and Sample_ID in manual_consensus_dict and gene in manual_consensus_dict[Sample_ID]:
				if args.verbose:
					print(f"Using manual consensus for Sample_ID {Sample_ID} gene {gene} from {manual_consensus_dict[Sample_ID][gene]}", file=sys.stderr)
			elif args.manual_consensus_file and args.verbose and Sample_ID in manual_consensus_dict:
				# Sample is in manual dict but this gene isn't
				print(f"Note: Sample_ID {Sample_ID} found in manual consensus file but gene {gene} not present. Available genes: {list(manual_consensus_dict[Sample_ID].keys())}", file=sys.stderr)
			
			if args.manual_consensus_file and Sample_ID in manual_consensus_dict and gene in manual_consensus_dict[Sample_ID]:
				try:
					consensus = read_in_fastaq_from_file(manual_consensus_dict[Sample_ID][gene])
					# Ensure the record.id matches the Sample_ID from metadata sheet
					consensus.id = Sample_ID
					consensus.description = Sample_ID + " manual_consensus"
					Clength = len(consensus)
					avqual = "NA"
					aliQual = "NA"
					nSeqs = "ManuallyCurratedConsensus"
					MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] = consensus
					MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] = Clength #read in the fasta/q sequence length from the file in manual_consensus_dict[Sample_ID][gene]
					MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] = avqual#round(sum(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"].letter_annotations["phred_quality"]) / len(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"]),1) if hasattr(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"], "letter_annotations") else "NA"
					MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["alignmentQual"] = aliQual
					MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["nSeqs"] = nSeqs
					write_consensus_fastq(MetaDict, args, Sample_ID, gene, consensus, consensus_fastq, aliQual, nSeqs)
					continue
				except Exception as e:
					print(f"Error: Failed to read manual consensus for Sample_ID {Sample_ID} gene {gene} from {manual_consensus_dict[Sample_ID][gene]}: {e}", file=sys.stderr)
					print(f"Falling back to ab1-based consensus generation for this sample.", file=sys.stderr)
					# Fall through to normal ab1 processing

			if args.resume and os.path.exists(consensus_fastq) and os.path.exists(consensus_fastq.replace(".fq", ".stats.txt")):
				if args.verbose:
					print(f"{consensus_fastq} already exists. Loading this file instead of re-creating it.", file=sys.stderr)
				consensus = SeqIO.read(consensus_fastq, "fastq")
				with open(consensus_fastq.replace(".fq", ".stats.txt"), "r") as CS:
					for line in CS:
						#print(line, file=sys.stderr)
						#input("Press Enter to continue...")  # Debugging pause
						aliQual, nSeqs = line.strip().split("\t")
						aliQual = float(aliQual)
						# Try to convert nSeqs to int, but keep as string if it's "ManuallyCurratedConsensus"
						try:
							nSeqs = int(nSeqs)
						except ValueError:
							pass  # Keep nSeqs as string (e.g., "ManuallyCurratedConsensus")
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] = consensus
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] = len(consensus) if consensus is not None else None
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] = round(sum(consensus.letter_annotations["phred_quality"]) / len(consensus),1) if consensus is not None else None
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["alignmentQual"] = aliQual
				MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["nSeqs"] = nSeqs
				#pprint.pprint(MetaDict[Sample_ID][gene]["Seqs"]["Cons"])
				#input("Press Enter to continue...")  # Debugging pause	
				continue
			#make the consensus fastq file from the ab1 files
			#check through the ab1 files and see how many pass qual and length filters.
			#If there are two or more, make a consensus sequence from the two best ones.
			#If there is one that passes, use that one as the consensus sequence.
			#If there are none that pass, then skip this sample and gene.
			#If there are more than two, then compare the qualities and lengths of the sequences and select the two best ones to make a consensus sequence.
			l = args.min_read_len #100 #minimum sequence length for alignmetn
			nseqs_in_cons = sum([1 for L in MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["seqLens"] if L >= l])
			#if number of lengths above 100 is >=2: choose best two (these have already been qual-trimmed)
			if nseqs_in_cons >= 2:
				#get the index of the two sequences with the highest average quality:
				indices = sorted(range(len(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["avQuals"])), key=lambda i: MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["avQuals"][i], reverse=True)[:2]
				#make a list of the SeqRecords for those indices:
				records_to_consensus = [MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["records"][i] for i in indices]
				consensus, aliQual = build_consensus_record(records_to_consensus, sample_id=Sample_ID, sample_description ="consensus")
				nSeqs = 2
				#pprint.pprint(consensus)
				#print(aliQual, nSeqs, sep="\t", file=sys.stderr)
				#input("Press Enter to continue...")  # Debugging pause
			#if there is only one long enough, use it. 	
			elif nseqs_in_cons == 1:
				consensus = MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["records"][MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["seqLens"].index(max(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["seqLens"]))]
				aliQual = 0
				nSeqs = 1
			#if number of seqLens above 100 is 0: skip this sample
			elif nseqs_in_cons == 0:
				if args.verbose:
					print(f"Warning: Sample_ID {Sample_ID} gene {gene} has no sequences longer than ", l, "bps. Skipping this sample and gene.", sep="", file=sys.stderr)
				consensus = None
				aliQual = 0
				nSeqs = 0
				
			#now put the consensus sequence into the metadata dictionary:
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] = consensus
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] = len(consensus) if consensus is not None else None 
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] = round(sum(consensus.letter_annotations["phred_quality"]) / len(consensus),1) if consensus is not None else None
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["alignmentQual"] = aliQual
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["nSeqs"] = nSeqs

			write_consensus_fastq(MetaDict, args, Sample_ID, gene, consensus, consensus_fastq, aliQual, nSeqs)


	print("	Samples processed: ", counter, end = "\r", file=sys.stderr, sep="")
	return MetaDict		

def write_consensus_fastq(MetaDict, args, Sample_ID, gene, consensus, consensus_fastq, aliQual, nSeqs):
				#write the consensus sequence to a fasta file
			if consensus is not None:
				consensus_fasta = MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"]
				if not os.path.exists(os.path.dirname(consensus_fasta)):
					os.makedirs(os.path.dirname(consensus_fasta), exist_ok=True)
				with open(consensus_fasta, "w") as fasta_handle:
					SeqIO.write(consensus, fasta_handle, "fasta")
				
				#write the consensus sequence to a fastq file
				if "phred_quality" in consensus.letter_annotations:
					if not os.path.exists(os.path.dirname(consensus_fastq)):
						os.makedirs(os.path.dirname(consensus_fastq), exist_ok=True)
					with open(consensus_fastq, "w") as fastq_handle:
						SeqIO.write(consensus, fastq_handle, "fastq")
				
				#write the consensus sequence stats to a file
				stats_file = consensus_fastq.replace(".fq", ".stats.txt")
				with open(stats_file, "w") as stats_handle:
					stats_handle.write(f"{aliQual}\t{nSeqs}\n")

def read_in_manual_consensus_file(args):
	print(f"Reading in manual consensus file: {args.manual_consensus_file}", file=sys.stderr)
	flag1 = True
	should_be_empty = ["", "", ""] #just to check that the columns that should be empty are actually empty.
	manual_consensus_dict = {}
	for line in open(args.manual_consensus_file, "r", encoding='utf-8-sig'):
		line = line.strip()
		# Skip empty lines
		if not line:
			continue
		if line.startswith("Sample_ID"):
			continue
		bits = line.split(",")
		# Skip malformed lines
		if len(bits) < 9:
			if args.verbose:
				print(f"Warning: Skipping malformed line in manual consensus file (expected 9 columns, got {len(bits)}): {line[:100]}", file=sys.stderr)
			continue
		sample_id = bits[0]
		# Skip lines with empty sample_id
		if not sample_id:
			continue
		plot = bits[1] #isn't needed in this file and can be ignored.
		directory = bits[2]
		ITSconsensus = bits[3]
		should_be_empty[0] = bits[4]
		ITS2consensus = bits[5]
		should_be_empty[1] = bits[6] 
		RBCLconsensus = bits[7]
		should_be_empty[2] = bits[8]
		if plot != "" and flag1:
			print(f"Notice: Manual consensus file {args.manual_consensus_file} has unexpected data in 'Plot' column. This is not necessary and will not be used here. Please assure that plots are correct in the standard metadata file as they will be used there.", file=sys.stderr)
			flag1 = False
		if should_be_empty != ["", "", ""]:
			print(f"Warning: Manual consensus file {args.manual_consensus_file} has unexpected data in reverse columns that should be empty: {should_be_empty}. Please only include consensus files in the forward primer positions. Exiting now.", file=sys.stderr)
			exit()
		if sample_id not in manual_consensus_dict and sum([ITS2consensus == "", RBCLconsensus == "", ITSconsensus == ""])<3: #if it is three, then don't bother adding in the sample ID -there are no consensuses for that sample.
			manual_consensus_dict[sample_id] = {}
			if ITSconsensus != "":
				manual_consensus_dict[sample_id]["ITS"] = args.prefix + directory + ITSconsensus
				if args.verbose:
					print(f"\tLoaded manual consensus for {sample_id} ITS: {args.prefix + directory + ITSconsensus}", file=sys.stderr)
			if ITS2consensus != "":
				manual_consensus_dict[sample_id]["ITS2"] = args.prefix + directory + ITS2consensus
				if args.verbose:
					print(f"\tLoaded manual consensus for {sample_id} ITS2: {args.prefix + directory + ITS2consensus}", file=sys.stderr)
			if RBCLconsensus != "":
				manual_consensus_dict[sample_id]["RBCL"] = args.prefix + directory + RBCLconsensus
				if args.verbose:
					print(f"\tLoaded manual consensus for {sample_id} RBCL: {args.prefix + directory + RBCLconsensus}", file=sys.stderr)
	print(f"\tLoaded manual consensus sequences for {len(manual_consensus_dict)} samples.", file=sys.stderr)
	#pprint.pprint(manual_consensus_dict, stream=sys.stderr)
	return manual_consensus_dict

def read_in_fastaq_from_file(filepath):
	"""
	Read in a fasta or fastq file and return a list of SeqRecords.
	"""
	if filepath.endswith(".fasta") or filepath.endswith(".fa") or filepath.endswith(".fna"):
		filetype = "fasta"
	elif filepath.endswith(".fastq") or filepath.endswith(".fq"):
		filetype = "fastq"
	records = []
	if not os.path.exists(filepath):
		print(f"Error: File {filepath} does not exist. Exiting now.", file=sys.stderr)
		exit()
	with open(filepath, "r") as handle:
		for record in SeqIO.parse(handle, filetype):
			records.append(record)
	if len(records) == 0:
		print(f"Error: No records found in file {filepath}.", file=sys.stderr)
		raise ValueError(f"Empty sequence file: {filepath}")
	if len(records) > 1:
		print(f"Warning: More than one record found in file {filepath}. Only the first record will be used.", file=sys.stderr)
	records = records[0]
	return records

def trim_seq_record(record, args, window_size=40):
	"""
	Trim a SeqRecord from both ends based on a moving average quality threshold.
	Returns a new trimmed SeqRecord, or a 5N/0Q record if too short after trimming.
	"""
	qual = args.qual if hasattr(args, 'qual') else 20  # Default quality threshold
	
	def mean(mylist): #tends to favour high values
		return sum(mylist)/len(mylist) if mylist else 0

	def geometric_mean(mylist): #tends to favour low values
		if not mylist:
			return 0
		product = 1
		for num in mylist:
			product *= num
		return product ** (1 / len(mylist))
	
	def harmonic_mean(mylist): #really favours low values
		if not mylist:
			return 0
		filtered = [x for x in mylist if x > 0]
		if not filtered:
			return 0
		return len(mylist) / sum(1 / x for x in mylist if x > 0)

	moving_average = 0
	i = 0
	j = 0
	length = len(record)

	#for old sanger sequencer machines, they didn't output the quality scores at all. So I need to skip this step if the quality scores are not present.
	if "phred_quality" not in record.letter_annotations:
		if args.verbose:
			# If the quality scores are not present, we cannot trim by quality.
			# We will return the original record without trimming.
			print(f"Warning: The sample {record.id} does not contain 'phred_quality' annotations. Cannot trim by quality.", file=sys.stderr)
		return record

	if sum(record.letter_annotations["phred_quality"]) == 0:
		if args.verbose:
			# If the quality scores are all zero, we cannot trim by quality.
			# We will return the original record without trimming.
			print(f"Warning: The sample {record.id} contains all zero 'phred_quality' annotations. Cannot trim by quality.", file=sys.stderr)
		return record

	if length > 5:
		while moving_average < qual and i < (length - window_size):
			moving_average = harmonic_mean(record.letter_annotations["phred_quality"][i:min(i+window_size, length)])
			i += 1
		moving_average = 0
		while moving_average < qual and j + i < (length - window_size):
			moving_average = harmonic_mean(record.letter_annotations["phred_quality"][max(0, length-j-window_size):length-j])
			j += 1
		if length - i - j > 5:
			trimmed = record[i:length-j]
			# --- Preserve abif_raw peak heights if present ---
			if "abif_raw" in record.annotations:
				trimmed.annotations = record.annotations.copy()
				# Optionally, you could trim PLOC2 to match the new sequence
				if "PLOC2" in trimmed.annotations["abif_raw"]:
					trimmed.annotations["abif_raw"]["PLOC2"] = record.annotations["abif_raw"]["PLOC2"][i:length-j]
			return trimmed
	# If too short or original length <= 5, return 5N/0Q record
	blank = record[:5]
	blank.seq = Seq("NNNNN")
	blank.letter_annotations = {"phred_quality": [0,0,0,0,0]}
	return blank

def trim_seq_by_peak_height(record, window_size=5, stddev_cutoff=0.2):
	"""
	Trim a SeqRecord from an ABI file based on average peak height.
	Bases are trimmed from both ends if the windowed average peak height is below (mean - stddev_cutoff * stddev).
	Returns a new trimmed SeqRecord.
	"""
	 # Create a blank record helper function
	def make_blank_record(rec):
		blank = SeqRecord(
			Seq("NNNNN"),
			id=rec.id,
			name=rec.name,
			description=rec.description
		)
		if "phred_quality" in rec.letter_annotations:
			blank.letter_annotations = {"phred_quality": [0,0,0,0,0]}
		return blank
	
	# Read ABI file
	#record = SeqIO.read(ab1_filepath, "abi")
	if record.seq is None or len(record.seq) < 6:
		return make_blank_record(record)
	
	if "abif_raw" not in record.annotations:
		print(f"Error: The record {record.id} does not contain 'abif_raw' annotations. Cannot process peak heights.", file=sys.stderr)
		return make_blank_record(record)
 	
	abif = record.annotations["abif_raw"]

	# Get peak heights for each base (A, C, G, T)
	a_peaks = np.array(abif["DATA9"], dtype=float)
	c_peaks = np.array(abif["DATA10"], dtype=float)
	g_peaks = np.array(abif["DATA11"], dtype=float)
	t_peaks = np.array(abif["DATA12"], dtype=float)

	# Get called base positions
	called_base_indices = abif["PLOC2"]
	seq = str(record.seq)
	peak_heights = []
	for i, base in enumerate(seq):
		idx = called_base_indices[i]
		if base == "A":
			peak_heights.append(a_peaks[idx])
		elif base == "C":
			peak_heights.append(c_peaks[idx])
		elif base == "G":
			peak_heights.append(g_peaks[idx])
		elif base == "T":
			peak_heights.append(t_peaks[idx])
		else:
			peak_heights.append(0.0)  # N or other
	#with open("./peaks.txt", "w") as out_fh:
		#for peak in peak_heights:
			#out_fh.write(f"{peak}\n")
	#pprint.pprint(peak_heights, stream=sys.stderr)
	#print(max(peak_heights), file=sys.stderr)
	#input("Press Enter to continue...")  # Debugging pause
	peak_heights = np.array(peak_heights)
	mean_peak = np.mean(peak_heights)
	std_peak = np.std(peak_heights)
	threshold = mean_peak - (stddev_cutoff * std_peak)
	#print(mean_peak, std_peak, threshold, sep="\t", file=sys.stderr)
	# Sliding window trim from both ends
	i, j = 0, 0
	length = len(peak_heights)
	# From left
	while i < (length - window_size):
		window_mean = np.mean(peak_heights[i:i+window_size])
		if window_mean >= threshold:
			break
		i += 1
	# From right
	while j < (length - window_size - i):
		window_mean = np.mean(peak_heights[length-j-window_size:length-j])
		if window_mean >= threshold:
			break
		j += 1

	# If too short after trimming, return a blank record
	if length - i - j < 5:
		return make_blank_record(record)

	trimmed = record[i:length-j]
	#print(f"Trimmed {i+j} bases leaving bases {i} to {length-j}", file=sys.stderr)
	#input("Press Enter to continue...")  # Debugging pause
	return trimmed

def build_consensus_record(rec_list, sample_id, sample_description="consensus"):
	# Align sequences 
	if len(rec_list) < 2:
		print("Error: Not enough sequences to build a consensus record. At least two sequences are required.", file=sys.stderr)
		raise ValueError("At least two sequences are required to build a consensus record.")
	elif len(rec_list) > 2:
		print("Warning: More than two sequences provided. Only the first two will be used for consensus building.", file=sys.stderr)
		rec_list = rec_list[:2]  # Use only the first two sequences

	rec1, rec2 = rec_list[0], rec_list[1]
	#print(rec1.seq, rec2.seq, sep="\n", file=sys.stderr)
	#("Press Enter to continue...")  # Debugging pause

	aligner = PairwiseAligner()

	# Set scoring similar to localxx: match=1, mismatch=0, no gap penalties
	aligner.mode = 'local'
	aligner.match_score = 1
	aligner.mismatch_score = 0
	aligner.open_gap_score = 0
	aligner.extend_gap_score = 0

	best_aln = aligner.align(rec1.seq, rec2.seq)[0] #might need to be rev1.seq
	#pprint.pprint(best_aln.format())
	#print("\n\n")
	if not best_aln:
		raise ValueError("No alignment found! Check %(sample_id)s for issues." % {"sample_id": sample_id})

	seq1_aln, seq2_aln = align_seqrecords_with_quality_from_alignment(best_aln, rec1, rec2)
	#print(rec1.seq)
	#print(rec2.seq)
	#print(seq1_aln.seq)
	#print(seq2_aln.seq)
	#input("Press Enter to continue...2")  # Debugging pause
	q1 = rec1.letter_annotations["phred_quality"]
	q2 = rec2.letter_annotations["phred_quality"]

	q1_aln = []
	q2_aln = []

	i, j = 0, 0  # indices for original quality arrays

	for b1, b2 in zip(seq1_aln, seq2_aln):
		if b1 != "-":
			q1_aln.append(q1[i])
			i += 1
		else:
			q1_aln.append(None)

		if b2 != "-":
			q2_aln.append(q2[j])
			j += 1
		else:
			q2_aln.append(None)

	consensus_seq = []
	consensus_qual = []

	for base1, base2, qval1, qval2 in zip(seq1_aln, seq2_aln, q1_aln, q2_aln):
		if base1 == base2 and base1 != "-":
			consensus_seq.append(base1)
			if qval1 is not None and qval2 is not None:
				consensus_qual.append(max(qval1, qval2))  # pick higher quality
			else:
				consensus_qual.append(qval1 or qval2 or 0)

		elif base1 == "-" and base2 != "-":
			consensus_seq.append(base2)
			consensus_qual.append(qval2 or 0)

		elif base2 == "-" and base1 != "-":
			consensus_seq.append(base1)
			consensus_qual.append(qval1 or 0)

		else:
			consensus_seq.append("N")
			consensus_qual.append(min(qval1 or 0, qval2 or 0))  # low confidence

	#pprint.pprint("".join(consensus_seq))
	#input("Press Enter to continue...3")  # Debugging pause
	return (SeqRecord(
		Seq("".join(consensus_seq)),
		id = sample_id,
		description = sample_description,
		letter_annotations = {"phred_quality": consensus_qual}
	), best_aln.score)

def align_seqrecords_with_quality_from_alignment(alignment, rec1, rec2):
	# Get aligned sequences with gaps
	seq1_aln = str(alignment[0])
	seq2_aln = str(alignment[1])

	q1 = rec1.letter_annotations.get("phred_quality", [])
	q2 = rec2.letter_annotations.get("phred_quality", [])

	q1_aln = []
	q2_aln = []

	i1 = i2 = 0
	for b1, b2 in zip(seq1_aln, seq2_aln):
		if b1 != "-":
			q1_aln.append(q1[i1] if i1 < len(q1) else 0)
			i1 += 1
		else:
			q1_aln.append(None)
		if b2 != "-":
			q2_aln.append(q2[i2] if i2 < len(q2) else 0)
			i2 += 1
		else:
			q2_aln.append(None)

	# Build aligned SeqRecords preserving metadata
	new_rec1 = SeqRecord(
		Seq(seq1_aln),
		id=rec1.id,
		name=rec1.name,
		description=rec1.description,
		letter_annotations={"phred_quality": [q if q is not None else 0 for q in q1_aln]},
		annotations=rec1.annotations,
		features=rec1.features,
	)
	new_rec2 = SeqRecord(
		Seq(seq2_aln),
		id=rec2.id,
		name=rec2.name,
		description=rec2.description,
		letter_annotations={"phred_quality": [q if q is not None else 0 for q in q2_aln]},
		annotations=rec2.annotations,
		features=rec2.features,
	)

	return new_rec1, new_rec2

#####CLASSIFY#######################################################
def classify_seqs_blastn(MetaDict, args):
	dbs = {"ITS": args.ITS_db, "ITS2": args.ITS2_db, "RBCL": args.RBCL_db}
	sls = {"ITS": args.species_list_ITS, "ITS2": args.species_list_ITS2, "RBCL": args.species_list_RBCL}

	gene_list = ["ITS", "ITS2", "RBCL"]
	if args.onlyITS:
		gene_list = ["ITS"]	
	
	for gene in gene_list:	
		print("\nClassifying sequences using BLASTN for ",  gene, ".",sep="",  file=sys.stderr)

		# Usage in your classify_seqs_blastn function:
		if sls[gene]:
			print("\tFiltering BLASTN database for species in the provided list.", file=sys.stderr)
			if not os.path.exists(sls[gene]):
				print(f"Error: Species list file {sls[gene]} not found. Re-run with the proper filename, or omit the species list. Exiting now.", file=sys.stderr)
				sys.exit(1)
			filtered_db = dbs[gene].replace(".fasta", "_filtered_") + os.path.basename(sls[gene].replace(".csv",".fasta"))
			#input(filtered_db)
			if os.path.exists(filtered_db) and args.resume: #CHANGE THIS BACK: NEEDS TO HAVE NO "NOT": including NOT is helping troubleshooting, but not needed in production
				print(f"\tFiltered BLASTN database {filtered_db} already exists. Using the previously filtered database.", file=sys.stderr)
			else:
				filter_fasta_by_species(dbs[gene], sls[gene], filtered_db)
			dbs[gene] = filtered_db #this replaces the db to be used with the filtered one

		#then make sure it has a blast database
		if not os.path.exists(dbs[gene] + ".nhr"):
			print("\tFormatting BLASTN database.", file=sys.stderr)
			#make the blastn database from the fasta file
			blastn_cmd = [
				"makeblastdb",
				"-in", dbs[gene],
				"-dbtype", "nucl",
			]
			try:
				subprocess.run(blastn_cmd, check=True)
			except subprocess.CalledProcessError as e:
				print(f"Error creating BLASTN database: {e}", file=sys.stderr)
				sys.exit(1)	

		#then on to using that db to classify the sequences. 
		counter = 0
		# Prepare jobs for parallel BLAST
		jobs = []
		for Sample_ID in MetaDict:
			#input(Sample_ID)
			consensus_record = MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"]
			if consensus_record is None or len(consensus_record) < 50:
				if args.verbose:
					print(f"Warning: Sample_ID {Sample_ID} gene {gene} has no consensus sequence or it is too short (<50 bps). Skipping this sample and gene.", file=sys.stderr)
				continue
			blastoutfile = os.path.join(args.output_dir, "02_classify_seqs", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_" + gene + ".blastnout")
			MetaDict[Sample_ID][gene]["Classifications"]["blastfile"] = blastoutfile
			query = MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"]
			if not os.path.exists(os.path.dirname(blastoutfile)):
				os.makedirs(os.path.dirname(blastoutfile), exist_ok=True)
			jobs.append((blastoutfile, MetaDict, Sample_ID, gene, dbs, args, query))

		# Run BLAST jobs in parallel
		with ThreadPoolExecutor(max_workers=args.cpus if args.cpus > 0 else 1) as executor:
			futures = [executor.submit(do_blast, *job) for job in jobs]
			for i, future in enumerate(as_completed(futures), 1):
				if i % 10 == 0:
					print(f"\tSamples processed: {i}", end="\r", file=sys.stderr)

	MetaDict = parse_blastn_results(MetaDict, gene_list, args)
	MetaDict = retry_blastn_on_trimmed_reads(MetaDict, gene_list, dbs, args)
	return MetaDict

def parse_blastn_results(MetaDict, gene_list, args):
	for gene in gene_list:
		print("\nParsing BLASTN output for ",  gene, ".",sep="",  file=sys.stderr)
		counter = 0
		for Sample_ID in MetaDict:			
			counter += 1
			if counter % 10 == 0:
				print("\tSamples processed: ", counter, end = "\r", file=sys.stderr, sep="")
			#read the blastn output file and parse it for the "best" hit
			blastfile = MetaDict[Sample_ID][gene]["Classifications"]["blastfile"]
			best_hit, consensus_taxonomy, classification_type = parse_blastn_output(MetaDict, Sample_ID, gene, args, blastfile) #this outputs the trimmed_list too, but we don't need it here
			#print("\t\tBest hit for ", Sample_ID, ":", best_hit, sep="", file=sys.stderr)
			SH_number = None
			if best_hit is not None and re.search("s__", consensus_taxonomy): #make sure that the best hit worked and goes to species level, if so, add in the SH number
				SH_number = best_hit.split("\t")[1].split("|")[2]  # Extract the SH number from the best hit
			MetaDict = update_MetaDict_with_blast_results(best_hit, consensus_taxonomy, MetaDict, Sample_ID, gene, args, blastfile, classification_type, SH_number)
		print("\tSamples processed: ", counter, end = "\r", file=sys.stderr, sep="")	

	#pprint.pprint(MetaDict["N230_11"]["ITS"])
	#input("Press Enter to continue...")  # Debugging pause
	#If there is no good blast hit, try again by blasting each trimmed read and see if either of those work. 
	return MetaDict

def retry_blastn_on_trimmed_reads(MetaDict, gene_list, dbs, args):	
	for gene in gene_list:
		print("\nRetrying the BLASTN for ",  gene, " for samples whose consensus had no hit. Using trimmed reads to see if they work better.",sep="",  file=sys.stderr)
		counter = 0
		for Sample_ID in MetaDict:
			#check if the SH_number is empty (i.e.: no species-level hit), if so, try to classify the trimmed reads
			if MetaDict[Sample_ID][gene]["Classifications"]["taxonomy"] == None:
				counter += 1
				if counter % 10 == 0:
					print("\tSamples processed: ", counter, end = "\r", file=sys.stderr, sep="")
				
				# need to write trimmed fqs to fa
				for fq in MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"]:
					#input(fq)
					if not os.path.exists(fq):
						#print("	not found")
						if args.verbose:
							print(f"Warning: Expected trimmed fastq file {fq} not found for Sample_ID {Sample_ID} gene {gene}. Skipping this sample and gene.", file=sys.stderr)
						continue
					SeqIO.convert(fq, "fastq", fq.replace(".fq", ".fa"), "fasta")  # Convert fastq to fasta for BLASTN

				# do the blastn
				blast_jobs = []
				blastfilelist = []
				for read, i in zip(["Forward", "Reverse"], range(len(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"]))):

					blastoutfile = os.path.join(args.output_dir, "02_classify_seqs", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_" + gene +"."+read+".blastnout")
					blastfilelist.append(blastoutfile)
					query = MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"][i].replace(".fq", ".fa")
					#print(query)
					#input(os.path.isfile(query))

					if not os.path.exists(query):
						#print("not found")
						if args.verbose:
							print(f"Warning: Expected read {query} not found for Sample_ID {Sample_ID} gene {gene}. Skipping this read.", file=sys.stderr)
						continue
					blast_jobs.append((blastoutfile, MetaDict, Sample_ID, gene, dbs, args, query))

				# Run BLAST jobs in parallel for this sample
				with ThreadPoolExecutor(max_workers=args.cpus if args.cpus > 0 else 1) as executor:
					futures = [executor.submit(do_blast, *job) for job in blast_jobs]
					for future in as_completed(futures):
						future.result()  # To raise exceptions if any
					
				# parse the blastn output as before	
				#blastfilelist = [MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"][j].replace("01_consensuses", "02_classify_seqs").replace(".fq", ".blastnout")for j in range(len(MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"]))]
				best_hits = []
				consensus_taxonomies = []

				for blastfile in blastfilelist:
					#need to go through those blast files and get the best hit from each, then compare the best hits, decide on which read should be used as the consensus and update everythinig.
					best_hit, consensus_taxonomy, classification_type = parse_blastn_output(MetaDict, Sample_ID, gene, args, blastfile) #
					best_hits.append(best_hit)
					consensus_taxonomies.append(consensus_taxonomy)

				#pprint.pprint(best_hits)
				#if both best hits are None, just continue without changing anything.
				if all(hit is None for hit in best_hits):
					#print("\t\tboth none", file=sys.stderr)
					continue  # No valid hits found, skip to next sample
				elif sum(hit is not None for hit in best_hits) == 1:
					#print("\t\t one not none", file=sys.stderr)
					#if exactly one of the hits is NOT none, use that one.
					index_of_best_hit = best_hits.index(next(hit for hit in best_hits))
				elif sum(hit is not None for hit in best_hits) >= 1:
					index_of_best_hit = max(range(len(best_hits)), key=lambda i: float(best_hits[i].split("\t")[2]) if best_hits[i] is not None else 0)
				
				best_hit = best_hits[index_of_best_hit]
				consensus_taxonomy = consensus_taxonomies[index_of_best_hit]
				blastfile = blastfilelist[index_of_best_hit]  # Update the blastfile to the one with the best hit
				#if Sample_ID =="N692_4":
				#	pprint.pprint(blastfilelist)
				#	input("how does blastfile list look for this sample?") 
				if best_hit is None:
					if args.verbose:
						print(f"Warning: Sample_ID {Sample_ID} gene {gene} has no valid BLASTN hits from trimmed reads. Skipping this sample and gene.", file=sys.stderr)
					continue
				SH_number = None
				if best_hit is not None and re.search("s__", consensus_taxonomy): #make sure that the best hit worked and goes to species level, if so, add in the SH number	
					SH_number = best_hit.split("\t")[1].split("|")[2]  # Extract the SH number from the best hit
				MetaDict = update_consensus_in_MetaDict(MetaDict, Sample_ID, gene, index_of_best_hit)
				MetaDict = update_MetaDict_with_blast_results(best_hit, consensus_taxonomy, MetaDict, Sample_ID, gene, args, blastfile, classification_type, SH_number)
				#pprint.pprint(MetaDict[Sample_ID][gene]["Classifications"]) 
				#pprint.pprint(MetaDict[Sample_ID][gene]["Seqs"]) 
				#input("Press Enter to continue...")  # Debugging pause	
		print("\tSamples processed: ", counter, end = "\r", file=sys.stderr, sep="")

	#pprint.pprint(MetaDict["N230_11"]["ITS"])
	#input("Press Enter to continue...")  # Debugging pause
	#	pprint.pprint(MetaDict["N692_4"]["ITS"])
	#	input("2! does it have a consensus fasta/q? YES quals")	
	
	return MetaDict

def update_consensus_in_MetaDict(MetaDict, Sample_ID, gene, index_of_best_hit):
	#if Sample_ID =="N692_4":
	#	pprint.pprint(MetaDict["N692_4"]["ITS"])
	#	input("3? does it have a consensus fasta/q?")	
	#update the MetaDict so the consensus reflects the best hit:
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] = MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["records"][index_of_best_hit]  # Update the consensus record with the best hit's record
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] = len(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"])
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] = round(sum(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"].letter_annotations["phred_quality"]) / len(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"]), 1) 
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["alignmentQual"] = 0  # Assuming no alignment quality for this case
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["nSeqs"] = 1
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fastq"] = MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"][index_of_best_hit]
	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"] = MetaDict[Sample_ID][gene]["Seqs"]["ab1s"]["fq_files"][index_of_best_hit].replace(".fq", ".fa")
	return MetaDict

def do_blast(blastoutfile, MetaDict, Sample_ID, gene, dbs, args, query):
	#MetaDict[Sample_ID][gene]["Classifications"]["blastfile"] = blastoutfile
	if not os.path.exists(os.path.dirname(blastoutfile)):
		os.makedirs(os.path.dirname(blastoutfile), exist_ok=True)
	if os.path.exists(blastoutfile) and args.resume:
		if args.verbose:
			print(f"Warning: BLASTN output file {blastoutfile} already exists. Skipping classification for Sample_ID {Sample_ID} gene {gene}.", file=sys.stderr)
		return 1  # Skip if the file already exists and skip_existing_blastn is True
	#make the blastn command
	blastn_cmd = [
		"blastn",
		"-query", query,
		"-db", dbs[gene],
		"-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand",
		"-max_target_seqs", "50",
		"-num_threads", "1", #str(args.cpus) if args.cpus > 0 else "1",
		"-out", blastoutfile
	]
	try:
		subprocess.run(blastn_cmd, check=True, capture_output=True)
	except subprocess.CalledProcessError as e:
		print(f"Error running BLASTN for Sample_ID {Sample_ID} gene {gene}: {e}", file=sys.stderr)

def filter_fasta_by_species(db_path, species_list_path, output_path):
	"""
	Filter a FASTA file by genus and species using a lookup dictionary from a species list.
	The species list should be a CSV/TSV with genus in the first column and species in the second.
	Only records whose genus and species are in the list are written to output_path.
	"""

	# Build the genus: set(species) dictionary
	species_dict = {}
	with open(species_list_path) as f:
		for line in f:
			if not line.strip():
				continue

			parts = line.strip().replace(",", "\t").replace(" ", "\t").replace("_", "\t").split("\t")  # handles CSV or TSV
			if len(parts) < 2:
				continue
			genus, species = parts[0], parts[1]
			genus = genus.lower()  # Normalize genus to lowercase
			species = species.lower()  # Normalize species to lowercase
			if genus not in species_dict:
				species_dict[genus] = set()
			species_dict[genus].add(species)
	#input("did spp dict get mad correctly?")
	#pprint.pprint(species_dict, stream=sys.stderr)  # Debugging line to see the species_dict
	# Filter the FASTA
	entries = 0
	with open(output_path, 'w') as out_fh:
		for record in SeqIO.parse(db_path, "fasta"):
			# Try to extract genus and species from record.id (assumes "Genus_species" or similar)
			# Adjust this parsing as needed for your FASTA header format!
			#input(record.id)
			taxonomy = record.id.split("|")[-1]
			#if re.search("rbcL", db_path):
			#	taxonomy = reformat_RBCL_taxonomy(taxonomy) #something to test if it's teh rbcl database
			taxa_levels = taxonomy.split(";")

			if len(taxa_levels) < 7:
				continue
			#print(record.id, file=sys.stderr)
			#print(taxonomy, file=sys.stderr)
			#print(taxa_levels, file=sys.stderr)
			genus = taxa_levels[5].lstrip("g__").lower() if len(taxa_levels) > 5 else None
			species = taxa_levels[6].split("_")[3].lower() if len(taxa_levels) > 6 else None
			#input(genus+" " + species)
			if genus in species_dict and species in species_dict[genus]:
				entries += 1
				SeqIO.write(record, out_fh, "fasta")
		print("\t\tFiltered FASTA database saved to", output_path, "with", entries, "entries matching the species list.", file=sys.stderr)

def parse_blastn_output(MetaDict, Sample_ID, gene, args, blastfile):
	#this needs to read in the blastn output file and parse it for the "best" hit
	#pprint.pprint(MetaDict[Sample_ID][gene]["Classifications"])
	#input("Press Enter to continue...")  # Debugging pause
	

	#print(blastfile_list, file=sys.stderr)
	if not blastfile:
		if args.verbose:
			print(f"Error: No BLASTN output file specified for Sample_ID {Sample_ID} gene {gene}.", file=sys.stderr)
		return None, None, None
	if not os.path.exists(blastfile):
		if args.verbose:
			print(f"Error: BLASTN output file {blastfile} does not exist for Sample_ID {Sample_ID} gene {gene}.", file=sys.stderr)
		return None, None, None
	
	with open(blastfile, "r") as blast_fh:
		#need to read in hits
		blast_list = [line.strip("\n") for line in blast_fh if int(line.split("\t")[3]) > args.minMatchLength and test_scov_qcov(line, args.minSCOV, args.minQCOV)]
	if len(blast_list) == 0:
		if args.verbose:
			print(f"Warning: No hits found for Sample_ID {Sample_ID} gene {gene} that pass the filters. Skipping this sample and gene.", file=sys.stderr)
		return None, None, None
	#if gene == "RBCL": #this is because the RBCL database is formatted oddly, I'm here changing some things around so that it's formatted correctly.
	#	blast_list = [reformat_RBCL_blast(hit) for hit in blast_list]

	if args.no_Incertae_Sedis:
		no_incertae_list = [hit for hit in blast_list if not re.search("sedis", hit)] #get rid of all hits with incertae sedis in the running for best hit
		if len(no_incertae_list) > 0: #if there's any left #
			blast_list = no_incertae_list 

	if args.adjustPidents:
		blast_list = [p_adjust(hit) for hit in blast_list]
		
	#find the best hit based on pident and length
	best_hit = max(blast_list, key=lambda x: (float(x.split("\t")[2]), int(x.split("\t")[3]))) #pick the highest pident BLAST hit breaking ties based on length
	#input("Press Enter to continue...")  # Debugging pause
	#find other hits that are within 5% of the best hit
	max_pident = float(best_hit.split("\t")[2])
	min_pident = max_pident - args.pidentDiff
	trimmed_list = [hit for hit in blast_list if float(hit.split("\t")[2]) >= min_pident]
	#if Sample_ID =="N230_11":
	#	print(f"Sample_ID {Sample_ID} gene {gene} has {len(trimmed_list)} hits after filtering.", file=sys.stderr)
	#	pprint.pprint(trimmed_list, stream=sys.stderr)  # Debugging line to see the blast hits
	#	input("Press Enter to continue...")  # Debugging pause

	if len(trimmed_list) == 0:
		return None, None, None  # No hits above the minimum pident threshold, return MetaDict unchanged
	
	consensus_taxonomy = get_consensus_taxonomy(trimmed_list, args)#returns the taxonomy string where all those hits agree. 
	#if Sample_ID =="N230_11":
	#	print(f"Sample_ID {Sample_ID} gene {gene} has {len(trimmed_list)} hits after filtering.", file=sys.stderr)
	#	pprint.pprint(consensus_taxonomy, stream=sys.stderr)  # Debugging line to see the blast hits
	#	input("Press Enter to continue...")  # Debugging pause

	#now need to truncate the taxonomy based on the pident of the best hit using args.pident_species etc.
	#based on the pident cutoffs, truncate consensus_taxonomy to the appropriate level
	pidents = [args.min_pident_species, args.min_pident_genus, args.min_pident_family, args.min_pident_order, args.min_pident_class, args.min_pident_phylum, args.min_pident_kingdom]
	consensus_taxonomy = truncate_taxonomy_by_pident(consensus_taxonomy, max_pident, pidents)
	
	if consensus_taxonomy is None:
		if args.verbose:
			print(f"Warning: No valid taxonomy found for Sample_ID {Sample_ID} gene {gene} with pident {max_pident}. Skipping this sample and gene.", file=sys.stderr)
		return None, None, None
	
	return(best_hit, consensus_taxonomy, "blastn")

def truncate_taxonomy_by_pident(taxonomy, pident, pident_thresholds):
	"""
	Truncate taxonomy string to the appropriate level based on pident and thresholds.
	taxonomy: string, e.g. "k__Fungi;p__Ascomycota;c__Sordariomycetes;..."
	pident: float, percent identity
	pident_thresholds: list of thresholds, most specific to least (species to kingdom)
	"""
	# Count how many thresholds are passed
	n_levels = sum(pident >= threshold for threshold in pident_thresholds)
	if n_levels > 0:
		return ";".join(taxonomy.split(";")[:n_levels])
	else:
		return None

def update_MetaDict_with_blast_results(best_hit, consensus_taxonomy, MetaDict, Sample_ID, gene, args, blastfile, classification_type, SH_number):
	#print(best_hit, file=sys.stderr)
	#print(MetaDict[Sample_ID][gene]["Classifications"]["best_hit"], file=sys.stderr)
	#get consensus taxonomy from the top hits in trimmed_list
	#need to look through all teh hits in trimmed list and see if they agree at each taxonomic level, report a taxonomy string truncated to where they agree.
	MetaDict[Sample_ID][gene]["Classifications"]["taxonomy"] = consensus_taxonomy
	#put best_hit into the MetaDict for this Sample_ID and gene
	MetaDict[Sample_ID][gene]["Classifications"]["best_hit"] = best_hit
	MetaDict[Sample_ID][gene]["Classifications"]["blastfile"] = blastfile
	MetaDict[Sample_ID][gene]["Classifications"]["classification_type"] = classification_type
	MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] = SH_number
	if gene =="ITS":
		MetaDict[Sample_ID][gene]["Clustering"]["cluster_data"]["cluster_id"] = SH_number
	
	#if classification_type is None:
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] = None  # No valid classification, set consensus record to None
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"] = None  # No valid classification, set consensus record to None
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fastq"] = None  # No valid classification, set consensus record to None
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] = 0  # No valid classification, set consensus record to None
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] = None  # No valid classification, set consensus record to None
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["alignmentQual"] = 0  # No valid classification, set consensus record to None
	#	MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["nSeqs"] = 0  # No valid classification, set consensus record to None

	return MetaDict

def test_for_pident(trimmed_list, pident_list, i):
	#needs to return True if the max pident in greater than the threshold of pident_list[i]	
	return(max([float(hit.split("	")[2]) for hit in trimmed_list]) > pident_list[i])

def test_scov_qcov(line, minSCOV, minQCOV):
	#make sure that scov>minscov and qcov>minqcov, but if one is close to 1, then it doesn't matter how low the other is:
	scov = calc_scov(line)
	qcov = calc_qcov(line)
	if abs(qcov - 1.0) < 0.05 or abs(scov - 1.0) < 0.05: #if one is within 5% of the matchlen, accept the blasthit.
		return(True)
	return(scov > minSCOV and qcov > minQCOV)

def calc_scov(hit):
	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,sstrand = hit.strip("\n").split("	")
	return(float(length) / float(slen))

def calc_qcov(hit):
	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,sstrand = hit.strip("\n").split("	")
	return(float(length) / float(qlen))

def reformat_RBCL_blast(hit):
	#should be defunct
	input("This should be defunct, but is still being called. Please fix this. 1")
	bits = hit.split("	")
	bits[1] = reformat_RBCL_taxonomy(bits[1])
	hit = "	".join(bits)
	return(hit)

def reformat_RBCL_taxonomy(bit):
	#should be defunct
	input("This should be defunct, but is still being called. Please fix this. 2")
	bit = bit.rstrip(";")
	taxonomy = bit.split(";")
	taxonomy = [re.sub(r'_[0-9]+$', '', tax) for tax in taxonomy]	
	spp = taxonomy[-1]
	spp = spp.lstrip("s__")
	bit = ";".join(taxonomy)
	bit = spp+"|"+spp+"|"+spp+"|reps|"+bit
	return(bit)

def get_consensus_taxonomy(trimmed_list, args):
	#this only needs to return a taxonomy string where the taxonomy at each level agrees (or agrees at >args.taxonomy_consensus_threshold).
	if args.pidentDiff == 0:
		# If pidentDiff is 0, we take the single hit with the highest evalue: I don't like this method as it's just taking the top blast hit without looking at any others. But it keeps being requested. I put a big warning in the help pages and made it not go this way by default - that's the best I can do. 
		maxEvalue = max(float(hit.split("\t")[10]) for hit in trimmed_list) 
		best_hit = next((hit for hit in trimmed_list if float(hit.split("\t")[10]) == maxEvalue), None)
		#print(best_hit)
		#print(best_hit.split("\t")[1].split("|")[-1] if best_hit else "None")
		return best_hit.split("\t")[1].split("|")[-1] if best_hit else "None"


	taxonomy = "None"
	#extract sseqid from blast hits, 
	sseqids = [hit.split("\t")[1] for hit in trimmed_list]
	#extract taxonomy from sseqid,
	taxonomy = [sseqid.split("|")[-1] for sseqid in sseqids]
	#pprint.pprint(taxonomy)
	#input("Press Enter to continue...")  # Debugging pause
	#then check for agreement at each taxonomic level.
	# Split each taxonomy string into a list
	split_taxa = [tax.split(";") for tax in taxonomy]
	# Zip together and compare each level
	consensus = []
	for taxa_tuple in zip(*split_taxa):
		if all(t == taxa_tuple[0] for t in taxa_tuple):
			consensus.append(taxa_tuple[0])
		else:
			break
	return ";".join(consensus)

def p_adjust(hit):
	qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, sstrand = hit.split("	")
	length = int(length)
	qstart = int(qstart)
	qend = int(qend)
	qlen = int(qlen)
	sstart = int(sstart)
	send = int(send)
	slen = int(slen)
	pident = float(pident)
	totlen = length
	if sstrand == "plus": #both strand in the positive orientation: compare starts with starts and ends with ends
		if sstart > 1 and qstart > 1: #an overhang on the front
			totlen +=  min(sstart, qstart)
		qdif = qlen - qend
		sdif = slen - send
		if qdif > 1 and sdif > 1: #an overhang on the back
			totlen += min(sdif, qdif)
	elif sstrand =="minus": #the subject is in the reverse orientation:
		#compare send with qstart
		sdif = slen - sstart #if sstrand is 'minus' then sstart > send. I.e: sstart = 400, send = 1
		if sdif > 1 and qstart > 1:
			totlen += min(sdif, qstart) 
		#compare sstart with qend
		qdif = qlen - qend
		if send > 1 and qdif > 1:
			totlen += min(send, qdif)
	else:
		print("sstrand is niether 'plus' nor 'minus': ", sstrand, file=sys.stderr)
		exit()
	
	adjpident = (pident * length) / totlen
	
# 	if and sstrand =="minus":
# 		print(slen, sstart, send, sstrand, sdif)
# 		print(slen, "---", sstart,"|||", send, "---1\n   ", sdif,"   ", sstart-send, "   ", send )
# 		print(qlen, qstart, qend, qdif)
# 		print(length, totlen)
# 		print(pident, adjpident)
# 		input("how does that look?")
	length = str(length)
	qstart = str(qstart)
	qend = str(qend)
	qlen = str(qlen)
	sstart = str(sstart)
	send = str(send)
	slen = str(slen)
	pident = str(pident)
	adjpident = str(adjpident)
	newhit = "	".join([qseqid, sseqid, adjpident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, sstrand])
	return(newhit)

#####CLUSTER#######################################################
def cluster_seqs_vsearch(MetaDict, args):
	#get needed seqs and make them into a fasta
	#go through metadict and look at the consensus quality and length, print all that pass the filters to a fasta file.
	gene_list = ["ITS"] #only ever cluster ITS2...   , "ITS2", "RBCL"]
	if args.cluster_plants_too:
			gene_list.append("ITS2")
			gene_list.append("RBCL")

	if args.onlyITS:
		gene_list = ["ITS"]

	#make a fasta file with all the sequences to cluster
	for gene in gene_list:
		if gene in ["ITS"]: #only do ITSx trimming if the ITS region is in there...
			print("\nFinding sequences to trim with ITSx for ", gene,".", sep="", file=sys.stderr)
			fasta_for_ITSx = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_seqs_for_ITSx.fasta")
			if not os.path.exists(os.path.dirname(fasta_for_ITSx)):
				os.makedirs(os.path.dirname(fasta_for_ITSx), exist_ok=True)
			with open(fasta_for_ITSx, "w") as itsx_handle:
				counter = 0
				for Sample_ID in MetaDict:
					#pprint.pprint(MetaDict[Sample_ID][gene])
					if not should_sample_be_queued_for_ITSx_trimming(MetaDict, Sample_ID, gene, args):
						continue
					counter += 1
					if counter % 10 == 0:
						print("\tSamples suitable for ITSx trimming: ", counter, end = "\r", file=sys.stderr, sep="")
					#add the sequence to the ITSx queue
					#if Sample_ID == "N718_9":
					#	pprint.pprint(MetaDict[Sample_ID][gene], stream=sys.stderr)  # Debugging line to see the MetaDict for this sample and gene
					#	input("Press Enter to continue...")  # Debugging pause
					MetaDict[Sample_ID][gene]["Clustering"]["ITSx_trim_queue"] = True
					# Ensure record.id matches Sample_ID before writing (important for manual consensus files)
					record_to_write = MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"]
					if record_to_write.id != Sample_ID:
						record_to_write.id = Sample_ID
					SeqIO.write(record_to_write, itsx_handle, "fasta")
				print("\tSamples suitable for ITSx trimming: ", counter, end = "\n", file=sys.stderr, sep="")
			if counter == 0:
				print(f"No sequences found that meet the criteria for ITSx trimming for gene {gene}. Skipping ITSx step.", file=sys.stderr)
				continue  # Skip to next gene if no sequences to trim
			outfile = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_ITSX_output")
			#now run ITSx on the sequences that are queued
			if not os.path.exists(fasta_for_ITSx):	
				print(f"Error: Fasta file for ITSx {fasta_for_ITSx} does not exist. Cannot proceed with ITSx.", file=sys.stderr)
				exit(1)
			
			ITSx_cmd = [
				"ITSx",
				"-i", fasta_for_ITSx,
				"-o", outfile,
				"-t", "A,F", #-t is the list of organism groups, see here: https://microbiology.se/publ/itsx_users_guide.pdf, "A"=plants, "F"=fungi
				"--allow_reorder", "F",
				"--graphical", "F",
				"--preserve", "T",
				#"--positions", "F",
				#"--only_full", "T",
				"--partial", "1",
				"--save_regions", "all",
				"--cpu", str(args.cpus) if args.cpus > 0 else "1",
			]
			#the output file we care about is called: ITS_ITSX_output.full_and_partial.fasta and is in os.path.join(args.output_dir, "03_cluster_seqs", gene+"_ITSX_output.full_and_partial.fasta")
			ITSx_output_fasta = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_ITSX_output.full_and_partial.fasta")
			if not os.path.exists(ITSx_output_fasta) or not args.resume:
				print("Running ITSx on ", counter, " sequences for ", gene, ".", sep="", file=sys.stderr)
				try:
					subprocess.run(ITSx_cmd, check=True, capture_output=True)
				except subprocess.CalledProcessError as e:
					print(f"Error running ITSx: {e}", file=sys.stderr)
			else:
				print("ITSx has already been run, using previous data for ", gene, ".", sep="", file=sys.stderr)
			
			#now parse the uc file and update the cluster information in the MetaDict
			MetaDict = load_ITSx_output(MetaDict, gene, args, ITSx_output_fasta)
			#input("ITSx trimming complete. Press Enter to continue...")  # Debugging pause
		elif gene in ["ITS2", "RBCL"]:
			#need to set the Send_to_Clustering flag for the RBCL samples similar to how it was set for ITSx above.
			for Sample_ID in MetaDict:
				MetaDict[Sample_ID][gene]["Clustering"]["Send_to_Clustering"] = should_sample_be_clustered(MetaDict, Sample_ID, gene, args)

		#now do the clustering with vsearch.
		#pprint.pprint(MetaDict["N718_9"][gene])
		#input("Press Enter to continue...PRECLUSTER")  # Debugging pause
	for gene in gene_list:
		num_seqs_to_cluster = sum([1 for Sample_ID in MetaDict if MetaDict[Sample_ID][gene]["Clustering"]["Send_to_Clustering"]])
		print("Clustering ", num_seqs_to_cluster, " sequences for ", gene,".", sep="", file=sys.stderr)
		#input("Press Enter to continue...")  # Debugging pause
		if not os.path.exists(os.path.join(args.output_dir, "03_cluster_seqs")):
			os.makedirs(os.path.join(args.output_dir, "03_cluster_seqs"), exist_ok=True)
		
		fasta_for_clustering = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_seqs_for_clustering.fasta")
		with open(fasta_for_clustering, "w") as cluster_handle:
			#thecount = 0
			for Sample_ID in MetaDict:
				#if teh ITSx results are good (present) for that sample, then add it to the clustering queue
				#if the ITSx results are bad (absent), mark it for salvaging.
				if not MetaDict[Sample_ID][gene]["Clustering"]["Send_to_Clustering"]:
					continue
				
				MetaDict[Sample_ID][gene]["Clustering"]["Clustered"] = True  # Mark this sample as clustered
				SeqIO.write(MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"], cluster_handle, "fasta")
			#print(thecount, " sequences written to ", fasta_for_clustering, file=sys.stderr)
		#now cluster those seqs using vsearch
		if not os.path.exists(fasta_for_clustering):	
			print(f"Error: Fasta file for clustering {fasta_for_clustering} does not exist. Cannot proceed with clustering.", file=sys.stderr)
			exit(1)
		
		outfile = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_denovoCentroids.fasta")
		UCcountfile = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_clusters.uc")
		
		vsearch_cmd = [
			"vsearch", 
			"--cluster_fast", fasta_for_clustering, 
			"--minseqlength", "90",
			"--id", str(args.percentIdentityForCluster/100),
			"--centroids", outfile, 
			"--threads", str(args.cpus) if args.cpus > 0 else "1", 
			"--strand", "both", 
			"--fasta_width", "0", 
			"--sizeout", 
			"--uc", UCcountfile
		]
		
		if not os.path.exists(outfile) or not args.resume:
			#print("\tRunning VSEARCH clustering for ", num_seqs_to_cluster, " sequences for ", gene, ".", sep="", file=sys.stderr)
			#input("Press Enter to continue...")  # Debugging pause
			try:
				subprocess.run(vsearch_cmd, check=True, capture_output=True)
			except subprocess.CalledProcessError as e:
				print(f"Error running VSEARCH: {e}", file=sys.stderr)
		else:
			print("\tVSEARCH clustering has already been run, using previous data for ", gene, ".", sep="", file=sys.stderr)
		#now parse the uc file and update the cluster information in the MetaDict
		MetaDict = load_and_assign_clusters(MetaDict, gene, args)
		print("\tdone", end="\n", file=sys.stderr)		
	#input("VSEARCH clustering complete. Press Enter to continue...")  # Debugging pause
	#pprint.pprint(MetaDict["N230_11"]["ITS"])
	#input("Press Enter to continue...")  # Debugging pause

	return MetaDict

def load_ITSx_output(MetaDict, gene, args, ITSx_output_fasta):
	#here we are reading in the ITSx output and updating the MetaDict with the results
	#the ITSx output is in the form of a fasta file with the sequences
	with open(ITSx_output_fasta, "r") as itsx_handle:
		#read in fasta record by record and extract seq records
		counter = 0
		for record in SeqIO.parse(itsx_handle, "fasta"):
			Sample_ID = record.id.split(" ")[0] #typically have the sampleid then a space and consensus or something else, only want the Sample_ID.
			ITSx_fasta_filepath = os.path.join(args.output_dir, "01_consensuses", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_" + gene + ".ITSxExtract.fasta")
			counter += 1
			record.id = Sample_ID + " " + "ITSx_extract"
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"] = ITSx_fasta_filepath
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fastq"] = None
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] = record
			MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] = len(record.seq)
			MetaDict[Sample_ID][gene]["Clustering"]["ITSx_trimmed"] = True
			MetaDict[Sample_ID][gene]["Clustering"]["Send_to_Clustering"] = True  # Mark this sample as ready for clustering
			#It'd be nice to write the ITSx fasta to 01_consesuses as well.
			if not os.path.exists(ITSx_fasta_filepath) or not args.resume: # #if the file doesn't exist and we are not resuming, write it out
				with open(ITSx_fasta_filepath, "w") as itsx_fh:
					SeqIO.write(record, itsx_fh, "fasta")
			print("\tITSx trimming complete for ", counter, " samples.", sep="", end = "\r", file=sys.stderr)
	print("\tITSx trimming successful for ", counter, " samples.", sep="", end = "\n", file=sys.stderr)
	return(MetaDict)

def should_sample_be_queued_for_ITSx_trimming(MetaDict, Sample_ID, gene, args):
	#check if the consensus sequence is classified to the species level and if it has a high enough quality
	# the output of this will go to ITSx for ITS extraction
	# do not cluster if: (i.e.: return false if)
		#blast succeeded to species level. 
		#there is no consensus sequence or it is too short or has low average quality

	if MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] is not None:
		return(False)  # Skip this sample if it has a valid SH number, already clustered: only cluster sequences that were not classified to the species level.
	
	if MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] is None or MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] is None or MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] is None:
		return(False) #these have no consensus or it's a bunch on N's - nothing to do for them. 
		
	if MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] < args.minSeqLengthForCluster:		#pprint.pprint(MetaDict[Sample_ID])
		return(False)  # Skip this sample if it doesn't meet the criteria

	return(True)

def should_sample_be_clustered(MetaDict, Sample_ID, gene, args):
	#check if the consensus sequence is long enough and has a high enough quality
	# the output of this will go to ITSx for ITS extraction
	# do not cluster if: (i.e.: return false if)
		#blast succeeded to species level. 
		#there is no consensus sequence or it is too short or has low average quality

	if MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] is not None:
		return(False)  # Skip this sample if it has a valid SH number, already clustered: only cluster sequences that were not classified to the species level.
	
	if MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] is None or MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] is None or MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["record"] is None:
		return(False) #these have no consensus or it's a bunch on N's - nothing to do for them. 
		
	if MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] < args.minSeqLengthForCluster or MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["avQual"] < args.minAverageQualityForCluster:		#pprint.pprint(MetaDict[Sample_ID])
		return(False)  # Skip this sample if it doesn't meet the criteria

	return(True)

def load_and_assign_clusters(MetaDict, gene, args):
	#define file with clusters in it
	cluster_filename = args.output_dir + "03_cluster_seqs/" + gene + "_clusters.uc"
	clustDict = {}
	#counter = counter2 = 0
	SH_tom_set = set([])
	#here we go through the cluster file and annotate every sample's cluster of the samples that were clustered. 
	with open(cluster_filename, "r") as clustFH:
		for line in clustFH:
			rtype,clustid,length,pident,strand,_,_,cigar,qseqid,sseqid = line.rstrip("\n").split("	")
			clustid = "denovo_"+clustid #make sure the cluster ID is prefixed with "denovo_" to avoid confusion with SH numbers
			#print(rtype, clustid, length, pident, strand, cigar, qseqid, sseqid, file=sys.stderr)
			#SH_tom_set.add("SH_tom_"+clustid)
			#qseqid = "_".join(qseqidlong.split("_")[:-2]) #these should be Sample_IDs. But come from what I've named each sequences back when I made fastas and consensuses from ab1s. Need to standardise fasta record names. 
			#sseqid = "_".join(sseqidlong.split("_")[:-2])
			#seqType = qseqidlong.split("_")[-1] #C for Consensus, also read1, read2, etc for resequenced samples and for lowQ samples that ended up using a read instead of the cons.
			if qseqid not in clustDict:
				clustDict[qseqid] = {"qseqid":qseqid, "clusterID":None, "centroid":None, "cigar":None, "pident":None, "length":None, "clustSize":None, "strand":None}

			if rtype == "H": #H is for a hit to the centroid
				#counter2 += 1
				clustDict[qseqid]["clusterID"] = clustid
				clustDict[qseqid]["length"] = int(length)
				clustDict[qseqid]["pident"] = float(pident)
				clustDict[qseqid]["strand"] = strand
				clustDict[qseqid]["cigar"] = cigar
				clustDict[qseqid]["centroid"] = sseqid #was sseqidlong....
				#if parse_cigar_string(cigar) >.2:
				#	counter += 1
			elif rtype == "C": #description of the centroid
				clustDict[qseqid]["clustSize"] = int(length) #the length in C hits is the number in the cluster
				clustDict[qseqid]["clusterID"] = clustid
				clustDict[qseqid]["centroid"] = qseqid #was long - shouldn't change from how it was set with H
				clustDict[qseqid]["pident"] = 100 #becuase is the centroid
				clustDict[qseqid]["strand"] = "+"
			elif rtype == "S": #centroid, redundant record but has the length of the centroid instead of the cluster size
				#counter2 += 1
				clustDict[qseqid]["length"] = int(length)
				clustDict[qseqid]["cigar"] = length+"M" #it's in a cluster alone - it's matches are 100%

	#pprint.pprint(clustDict)
	#This bit adds in the cluster size for the samples which weren't centroids
	for SampleID in clustDict:
		#pprint.pprint(clustDict[SampleID])
		if clustDict[SampleID]["clustSize"] is None:
			centroid = clustDict[SampleID]["centroid"]
			clustDict[SampleID]["clustSize"] = clustDict[centroid]["clustSize"]
			
			
		#get the centroid data moved into position
		MetaDict[SampleID][gene]["Clustering"]["cluster_data"] = clustDict[SampleID]
		if MetaDict[SampleID][gene]["Classifications"]["SH_number"] is None:  
			MetaDict[SampleID][gene]["Classifications"]["SH_number"] = clustDict[SampleID]["clusterID"] 
		else:
			pprint.pprint(MetaDict[SampleID][gene])
			print("for some reason the SH_number is not None for ", SampleID, " in gene ", gene, ". This should not happen as samples with SH numbers should not have been passed to clustering. Exiting Now, Tom needs to fix this.", file=sys.stderr)
			exit(1)
	return(MetaDict)

#####SALVAGE#######################################################
def salvage_lowQ_seqs(MetaDict, args):
	#this will go through the MetaDict and find the samples that have attempedSalvaging set to True, and then try to salvage them by blasting a sequence of theirs to the clustering centroids
	#for each sample that has True under attemptedSalvaging, it will blast the consensus sequence to the centroids and see if it can find a match.
	gene_list = ["ITS"] #only ever salvage ITS for fungi....   , "ITS2", "RBCL"]
	if args.onlyITS:
		gene_list = ["ITS"]

	#need to make blastdb out of the clustered sequences:
	
	for gene in gene_list:
		centroids_fasta = os.path.join(args.output_dir, "03_cluster_seqs", gene+"_denovoCentroids_and_SHs.fasta")
		if not os.path.exists(centroids_fasta) or not args.resume:
			print(f"\tCentroids fasta file {centroids_fasta} does not exist. Creating it now with the centroids and the SH hits.", file=sys.stderr)
			#get list of SH numbers from the MetaDict
			SH_numbers = set()
			for Sample_ID in MetaDict:
				if MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] is not None:
					SH_numbers.add(MetaDict[Sample_ID][gene]["Classifications"]["SH_number"])
			#extract those sequences from UNITE and put them into a fasta file
			#print("SH_numbers:", SH_numbers, file=sys.stderr)
			#input("Press Enter to continue...")  # Debugging pause
			with open(centroids_fasta, "w") as centroids_handle:
				#write the denovo centroids first
				with open(os.path.join(args.output_dir, "03_cluster_seqs", gene+"_denovoCentroids.fasta"), "r") as denovo_handle:
					for record in SeqIO.parse(denovo_handle, "fasta"):
						SeqIO.write(record, centroids_handle, "fasta")
				#then write the SH numbers
				with open(args.ITS_db, "r") as unite_handle:
					for record in SeqIO.parse(unite_handle, "fasta"):
						#check if the record id starts with any of the SH numbers
						SH_number = record.id.split("|")[2]
						#print(SH_number)
						#input("correct?")  # Get the SH number from the record ID
						if SH_number in SH_numbers:
							#write the record to the centroids fasta file
							SeqIO.write(record, centroids_handle, "fasta")

			#put UNITE sequences into a fasta file along wih the denovo centroids.
		
		#make a blast db out of the centroids fasta file
		print("Creating BLASTN database from clustered sequences for salvaging low quality sequences.", file=sys.stderr)
		try:
			subprocess.run(["makeblastdb", "-in", centroids_fasta, "-dbtype", "nucl", "-out", centroids_fasta], check=True, capture_output=True)
		except subprocess.CalledProcessError as e:
			print(f"Error creating BLASTN database: {e}", file=sys.stderr)
			exit(1)
		print("\tdone", file=sys.stderr)
		counter = 0
		
		print("Salvaging low quality sequences for ",  gene, ".",sep="",  file=sys.stderr)
		for Sample_ID in MetaDict:
			if not should_sample_be_salvaged(MetaDict, Sample_ID, gene, args):
				continue
			counter += 1
			print("\tProcessing files:", counter, end="\r", file=sys.stderr)
			MetaDict[Sample_ID][gene]["Salvaging"]["attempted"] = True  # Mark this sample as not needing ITSx trimming anymore
			#pprint.pprint(MetaDict[Sample_ID][gene]["Seqs"])
			#input("Press Enter to continue...")  # Debugging pause
			do_salvage_blast(
				blastoutfile=os.path.join(args.output_dir, "02_classify_seqs", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_" + gene + ".salvage.blastnout"),
				Sample_ID=Sample_ID,
				gene=gene,
				dbs={"ITS":centroids_fasta, "ITS2":centroids_fasta, "RBCL":centroids_fasta},
				args=args,
				query=MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["fasta"]
			)

			#print("\n\n##################################################", file=sys.stderr)
			#print(Sample_ID)
			#print(MetaDict[SampleID]["Metadata"]["SalvageGroup"], file=sys.stderr)
			#print(MetaDict[SampleID][gene]["Classifications"]["taxonomy"], file=sys.stderr)

			OTU = parse_salvage_blastn_output(MetaDict, Sample_ID, gene, args, os.path.join(args.output_dir, "02_classify_seqs", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_" + gene + ".salvage.blastnout"))
			MetaDict = update_MetaDict_with_salvage_blast_results(MetaDict, Sample_ID, gene, OTU, args)
			#print out the Sample_ID and taxonomy, and whether clustering was attempted, and whether it was salvaged.
			#print("Centroid is", cluster_centroid, file=sys.stderr)
			#if cluster_centroid is not None:
			#	pprint.pprint(MetaDict[cluster_centroid], stream=sys.stderr)
			#input("Salvaging complete for sample. Press Enter to continue...")  # Debugging pause
		print("\tProcessing files:", counter, end="", file=sys.stderr)

	return(MetaDict)

def should_sample_be_salvaged(MetaDict, Sample_ID, gene, args):
	
	#MetaDict[Sample_ID][gene]["Clustering"]["ITSx_trim_queue"]
	#two cases when a sample should be salvaged:
	#1. The sample has >min_read_len_salvage bps in it's consensus, but did not class to species level and was not passed to ITSx for trimming.
	passed1 = 0
	if MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] is None: #must not be to species level
		passed1 += 1 
	if MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] is not None and MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["seqLen"] >= args.min_read_len_salvage: #must have a consensus sequence that is at least 150bps long
		passed1 += 1
	if MetaDict[Sample_ID][gene]["Clustering"]["ITSx_trim_queue"] is False: #must not have been passed to ITSx for trimming
		passed1 += 1
	
	if passed1 == 3: #if all three conditions are met, then the sample should be salvaged
		return(True)

	#2. The sample has ITSx_trim_queue set to True, meaning it was poor enough to be sent to ITSx, but has ITSx_trimmed set to False, meaning it didn't pass the ITSx trimming.
	if MetaDict[Sample_ID][gene]["Clustering"]["ITSx_trim_queue"] is True and MetaDict[Sample_ID][gene]["Clustering"]["ITSx_trimmed"] is False:
		return(True)
	
	return(False)

def do_salvage_blast(blastoutfile, Sample_ID, gene, dbs, args, query):
	#MetaDict[Sample_ID][gene]["Classifications"]["blastfile"] = blastoutfile
	if not os.path.exists(os.path.dirname(blastoutfile)):
		os.makedirs(os.path.dirname(blastoutfile), exist_ok=True)
	#not this bit - in case the clusters have changed, we need to re-run the blastn.
	#if os.path.exists(blastoutfile) and args.resume:
	#	if args.verbose:
	#		print(f"Warning: BLASTN output file {blastoutfile} already exists. Skipping classification for Sample_ID {Sample_ID} gene {gene}.", file=sys.stderr)
	#	return 1  # Skip if the file already exists and skip_existing_blastn is True
	#make the blastn command
	blastn_cmd = [
		"blastn",
		"-query", query,
		"-db", dbs[gene],
		"-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand",
		"-max_target_seqs", "50",
		"-num_threads", str(args.cpus) if args.cpus > 0 else "1",
		"-out", blastoutfile
	]
	try:
		subprocess.run(blastn_cmd, check=True, capture_output=True)
	except subprocess.CalledProcessError as e:
		print(f"Error running BLASTN for Sample_ID {Sample_ID} gene {gene}: {e}", file=sys.stderr)

def parse_salvage_blastn_output(MetaDict, Sample_ID, gene, args, blastfile):
	#this needs to read in the blastn output file and parse it for the "best" hit
	#pprint.pprint(MetaDict[Sample_ID][gene]["Classifications"])
	#input("Press Enter to continue...")  # Debugging pause
	

	#print(blastfile_list, file=sys.stderr)
	if not blastfile:
		if args.verbose:
			print(f"Error: No BLASTN output file specified for Sample_ID {Sample_ID} gene {gene}.", file=sys.stderr)
		return None
	if not os.path.exists(blastfile):
		if args.verbose:
			print(f"Error: BLASTN output file {blastfile} does not exist for Sample_ID {Sample_ID} gene {gene}.", file=sys.stderr)
		return None
	with open(blastfile, "r") as blast_fh:
		#need to read in hits
		blast_list = [line.strip("\n") for line in blast_fh if test_scov_qcov(line, args.minSCOV, args.minQCOV)] #only going to test qcov/scov here
	
	if len(blast_list) == 0:
		if args.verbose:
			print(f"Warning: No hits found for Sample_ID {Sample_ID} gene {gene} that pass the filters. Skipping this sample and gene.", file=sys.stderr)
		return None


	#find the best hit based on pident and length
	#best_hit = max(blast_list, key=lambda x: (float(x.split("\t")[11]))) #pick the BLAST hit with the highest eval
	# Find the maximum value for the key
	max_eval = max(float(x.split("\t")[11]) for x in blast_list)
	# Get all entries with that value
	best_hits = [x for x in blast_list if float(x.split("\t")[11]) == max_eval]
	
	#print("Best hit for Sample_ID", Sample_ID, ":", file=sys.stderr)
	#pprint.pprint(best_hits)	
	#input("Press Enter to continue...")  # Debugging pause
	#find other hits that are within 5% of the best hit
	cluster_centroids = [best_hit.split("\t")[1].split(";")[0] for best_hit in best_hits]  # The sseqid of the best hit is the cluster centroid
	#pprint.pprint(cluster_centroids, stream=sys.stderr)  # Debugging line to see the cluster centroids
	#input("Press Enter to continue...")  # Debugging pause
	cluster_centroids = [centroid.split("|")[2] if "|" in centroid else centroid for centroid in cluster_centroids ]  # Extract the SH number from the sseqid if it exists, otherwise use the whole string: becuase centroids are named by teh sampleid of the centorid whereas the SH ones are named with the OTU name. Will eventually need to translate the CCids to denovo_* ids and can then treat them similarly.
	#convert sampleIDs in the cluster_centroids to OTU names (i.e.: denovo_* which are treated like SH numbers)
	#pprint.pprint(cluster_centroids, stream=sys.stderr)  # Debugging line to see the cluster centroids
	
	#this is really what these are: OTUs. cluster centroids was a holdover from when we salvaged only to the clustering centroids.
	OTUs = [centroid if centroid.startswith("SH") else MetaDict[centroid][gene]["Clustering"]["cluster_data"]["clusterID"] for centroid in cluster_centroids]

	chosen_OTU = check_salvage_group(MetaDict, Sample_ID, gene, OTUs) #this is both checking if the cc is in the same sg, AND allowing for a the first of a couple similar hits to match. 
	if chosen_OTU is not None:
		return chosen_OTU
	return(None)
	
def check_salvage_group(MetaDict, focal_Sample_ID, gene, OTUs):
	"""
	For each cluster_centroid in the list, check if any member of that cluster
	is in the same salvage group as the focal_Sample_ID.
	Return the first such cluster_centroid found, or None if none found.
	"""
	#get the SH/CC of the focal sample (will be another sample as the cluster centroid, or an SH number depending on the blast hits) these need to be treated differently: 
	#the SH needs a lookup of all samples that share that SH and then a lookup of those samples' SG. 
	#The CC needs a lookup of the ClustID (denovo_*) and then that can be treated like the SH above.
	#scan through MetaDict for all samples in the same salvage group as the focal sample

	focal_SG = MetaDict[focal_Sample_ID]["Metadata"]["SalvageGroup"]
	#now for each OTU, extract all SGs that match that OTU:
	for OTU in OTUs:
		SGs = set()
		for Sample_ID in MetaDict:
			if MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] == OTU:
				SGs.add(MetaDict[Sample_ID]["Metadata"]["SalvageGroup"])
		if focal_SG in SGs:
			return OTU
	return None

def update_MetaDict_with_salvage_blast_results(MetaDict, Sample_ID, gene, OTU, args):
	#this updates the MetaDict with the results of the salvage blast
	if OTU is None:
		if args.verbose:
			print(f"Warning: No valid OTU found for Sample_ID {Sample_ID} gene {gene}. Skipping this sample and gene.", file=sys.stderr)
		return MetaDict
	#pprint.pprint(MetaDict[Sample_ID][gene]["Classifications"])
	#input("Press Enter to continue...")  # Debugging pause

	#need to update metadict now, these first are updated regardless of whether the OTU is an SH number or a denovo cluster.

	MetaDict[Sample_ID][gene]["Classifications"]["SH_number"] = OTU  # Set the SH number to the OTU found
	MetaDict[Sample_ID][gene]["Classifications"]["blastfile"] = os.path.join(args.output_dir, "02_classify_seqs", add_in_dirDepth(Sample_ID, gene, args), Sample_ID + "_" + gene + ".salvage.blastnout")
	MetaDict[Sample_ID][gene]["Classifications"]["classification_type"] = "salvage_blastn"
	MetaDict[Sample_ID][gene]["Clustering"]["Clustered"] = True  # Mark this sample as clustered
	MetaDict[Sample_ID][gene]["Clustering"]["clustered_by_salvage"] = True  # Mark this sample as clustered by salvaging
	MetaDict[Sample_ID][gene]["Salvaging"]["Successfull"] = True  #

	#if the OTU is an SH number, we can give it the whole taxonomy of that SH number.
	centroid = None

	if OTU.startswith("SH"):
		MetaDict[Sample_ID][gene]["Classifications"]["taxonomy"] = get_taxonomy_from_SH(MetaDict, OTU)
		MetaDict[Sample_ID][gene]["Classifications"]["best_hit"] = get_best_hit_from_SH(MetaDict, Sample_ID, gene, OTU) 
		
	#if the OTU is a denovo, we need to look up the centroid of that cluster and use it's taxonomy to classify the sample. 
	elif OTU.startswith("denovo_"):
		for tmp in MetaDict:
			if tmp == Sample_ID:  # Skip the sample itself
				continue
			if MetaDict[tmp][gene]["Classifications"]["SH_number"] == OTU:
				centroid = tmp
				break
	if centroid is not None:
		# If we found a centroid, we can use it
		MetaDict[Sample_ID][gene]["Classifications"]["taxonomy"] = MetaDict[centroid][gene]["Classifications"]["taxonomy"]
		MetaDict[Sample_ID][gene]["Classifications"]["best_hit"] = get_best_hit_from_SH(MetaDict, Sample_ID, gene, centroid) #yes, this needs to be centroid as thats what is in the blast file sseqid. Before it needed to be an SH number. 
		MetaDict[Sample_ID][gene]["Clustering"]["cluster_data"] = MetaDict[centroid][gene]["Clustering"]["cluster_data"]  # Use the cluster data from the centroid
	
	return MetaDict

def get_taxonomy_from_SH(MetaDict, OTU):
	#this will return the taxonomy of the SH number from the MetaDict
	#all samples with that SH will have the same taxonomy, so we can just return the first one we find.
	for Sample_ID in MetaDict:
		if MetaDict[Sample_ID]["ITS"]["Classifications"]["SH_number"] == OTU:
			return MetaDict[Sample_ID]["ITS"]["Classifications"]["taxonomy"]
	return None  # If no sample found with that SH number, return None

def get_best_hit_from_SH(MetaDict, Sample_ID, gene, OTU):
	#this will return the best hit of the SH number from the MetaDict
	#needs to find the blast file and parse it for the top hit that matches the OTU, then return that hit.
	with open(MetaDict[Sample_ID][gene]["Classifications"]["blastfile"], "r") as blast_fh:
		for line in blast_fh:
			if re.search(OTU, line):
				return line.strip("\n")  # Return the entire hit.
	return None

#####FUNGUILD#######################################################
def run_funguild(MetaDict, args):
	#runs the funguild pipeline - first needs to make the input file for the script
	gene = "ITS" #funguild is only for fungi

	#first make the input file: one column is sampleID, the next is the taxonomy
	input_file = write_FUNGUILD_input_file(MetaDict, gene, args)
	FUNGuild_executable = "Guilds_v1.1.py"  # Default assumes it's in PATH
	if args.FUNGuild_executable != "":
		#check if filepath exists
		if not os.path.exists(args.FUNGuild_executable):
			print("The path to the FUNGuild executable that you supplied is incorrect, I can find no file at:", args.FUNGuild_executable, file=sys.stderr)
			exit(1)
		FUNGuild_executable = args.FUNGuild_executable
	
#	funguild_cmd = [FUNGuild_executable, "-otu", input_file]
	funguild_cmd = [sys.executable, FUNGuild_executable, "-otu", input_file]
	try:
		subprocess.run(funguild_cmd, check=True, capture_output=True)
	except subprocess.CalledProcessError as e:
		print(f"Error running FUNGuild: {e}", file=sys.stderr)
		print("STDOUT:", e.stdout.decode(), file=sys.stderr)
		print("STDERR:", e.stderr.decode(), file=sys.stderr)
		exit()

	MetaDict = load_FUNGuild(MetaDict, gene, input_file)
	
	return(MetaDict)	

def write_FUNGUILD_input_file(MetaDict, gene, args):
	FUNGuildInputFile = args.output_dir + "04_FUNGuild/" + gene + "_FUNGuild.csv"
	with open (FUNGuildInputFile, "w") as FGFH:
		print("SampleID,taxonomy", file=FGFH)
		for SampleID in MetaDict:
				if MetaDict[SampleID][gene]["Classifications"]["taxonomy"] != "NoGoodHit" and MetaDict[SampleID][gene]["Classifications"]["taxonomy"] != "" and MetaDict[SampleID][gene]["Classifications"]["taxonomy"] != "NA":
					print(SampleID, MetaDict[SampleID][gene]["Classifications"]["taxonomy"], sep=",", file=FGFH)
	return(FUNGuildInputFile)			

def load_FUNGuild(MetaDict, gene, input_file):
	FUNGuildOutputFile = input_file.replace(".csv", ".guilds.txt")
	#print(input_file, "\n", FUNGuildOutputFile)
	with open(FUNGuildOutputFile, "r") as FGFH:
		lc = 0
		for line in FGFH:
			lc += 1
			if lc >1:
				line=line.replace("-", "None")  # replace dashes with None
				SampleID,taxonomy,Taxon,Taxon_Level,Trophic_Mode,Guild,Growth_Morphology,Trait,Confidence_Ranking,Notes,Citation_Source = line.rstrip("\n").split("	")
				MetaDict[SampleID][gene]["FUNGuild"] = {"Taxon":Taxon,"Taxon_Level":Taxon_Level,"Trophic_Mode":Trophic_Mode,"Guild":Guild,"Growth_Morphology":Growth_Morphology,"Trait":Trait,"Confidence_Ranking":Confidence_Ranking,"Notes":Notes,"Citation_Source":Citation_Source}
	return(MetaDict)

#####OUTPUTS#######################################################
def print_final_summary_Stats(MetaDict, args, FH, command_run):
	#this is just for some on-screen summary statistics at the end of the call, 
	#it should all be recoverable from the Sanger_Summary.txt file too. 
	#pprint.pprint(MetaDict)

	version = "unknown"
	version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "VERSION")
	with open(version_file) as f:
		version = f.read().strip()
	print("#", command_run, " # version: ", version, end="\n", sep="", file=FH)

	
	gene_list = ["ITS", "ITS2", "RBCL"]
	if args.onlyITS:
		gene_list = ["ITS"]

	db_used = {
		"ITS":os.path.basename(args.ITS_db),
		"ITS2":os.path.basename(args.ITS2_db),
		"RBCL":os.path.basename(args.RBCL_db)
	}

	if args.species_list_ITS:
		db_used["ITS"] = db_used["ITS"].replace(".fasta", "_"+os.path.basename(args.species_list_ITS).replace(".csv", ".fasta"))
	if args.species_list_ITS2:
		db_used["ITS2"] = db_used["ITS2"].replace(".fasta", "_"+os.path.basename(args.species_list_ITS2).replace(".csv", ".fasta"))
	if args.species_list_RBCL:
		db_used["RBCL"] = db_used["RBCL"].replace(".fasta", "_"+os.path.basename(args.species_list_RBCL).replace(".csv", ".fasta"))

	for gene in gene_list:
		#cluster_id = []
		classTable = {"Species":0, "Genus":0, "Family":0, "Order":0, "Class":0, "Phylum":0, "Kingdom":0, "unclassified":0}
		clustTable = {}
		clustTableHist = {}
		OTUTable = {}
		OTUTableHist = {}
		fungTable = {"Highly Probable":0, "Probable":0, "Possible":0, "-":0}
		fungTableGuild = {}
		funguildAttempts = 0
		for SampleID in MetaDict:

			if MetaDict[SampleID][gene]["Classifications"]["taxonomy"] is None:
				classTable["unclassified"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 7:
				classTable["Species"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 6:
				classTable["Genus"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 5:
				classTable["Family"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 4:
				classTable["Order"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 3:
				classTable["Class"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 2:
				classTable["Phylum"] += 1
			elif len(MetaDict[SampleID][gene]["Classifications"]["taxonomy"].split(";")) == 1:
				classTable["Kingdom"] += 1
			else:
				classTable["unclassified"] += 1

			
			if gene == "ITS":
				confRank = MetaDict[SampleID][gene]["FUNGuild"]["Confidence_Ranking"]
				guild = MetaDict[SampleID][gene]["FUNGuild"]["Guild"]
				if confRank != "NA" and confRank != "None":
					funguildAttempts += 1 
					if confRank not in fungTable:
						fungTable[confRank] = 0
					fungTable[confRank] += 1
				if guild not in ["-", "NA", "None"]:
					if guild not in fungTableGuild:
						fungTableGuild[guild] = 0
					fungTableGuild[guild] += 1

			SH_id = MetaDict[SampleID][gene]["Classifications"]["SH_number"] 		
			if gene == "ITS" or args.cluster_plants_too:			
				cluster_id = MetaDict[SampleID][gene]["Clustering"]["cluster_data"]["clusterID"]
				clustered = MetaDict[SampleID][gene]["Clustering"]["Clustered"]
				#pprint.pprint(MetaDict[SampleID][gene]["Clustering"])
				#input("Press Enter to continue...")  # Debugging pause	
				#print(cluster_id, file=sys.stderr)
				#input("Press Enter to continue...")  # Debugging pause
				#clustered samples are those that have a clustered as true and clustered_by_salvage as false:
				if clustered == True and MetaDict[SampleID][gene]["Clustering"]["clustered_by_salvage"] == False:
					if cluster_id not in clustTable:
						clustTable[cluster_id] = 0
					clustTable[cluster_id] += 1

			if SH_id is not None:
				if SH_id not in OTUTable:
					OTUTable[SH_id] = 0
				OTUTable[SH_id] += 1
		if gene == "ITS" or args.cluster_plants_too:
			for cluster_id in clustTable:
				size = clustTable[cluster_id]
				if size not in clustTableHist:
					clustTableHist[size] = 0
				clustTableHist[size] +=1
	
		for SH_id in OTUTable:
			size = OTUTable[SH_id]
			if size not in OTUTableHist:
				OTUTableHist[size] = 0
			OTUTableHist[size] +=1
		
		minClusterSize = "NA"
		maxClusterSize = "NA"
		aveClusterSize = "NA"
		medClusterSize = "NA"
		numSingletonClusts = "NA"
		
		minOTUClusterSize = "NA"
		maxOTUClusterSize = "NA"
		aveOTUClusterSize = "NA"
		medOTUClusterSize = "NA"
		numOTUSingletonClusts = "NA"
			
		if len(clustTableHist) > 0:
			
			minClusterSize = min(clustTableHist.keys())
			maxClusterSize = max(clustTableHist.keys())
			aveClusterSize = round(sum(clustTable.values())/len(clustTable),2)
			medClusterSize = median(clustTable.values())
			numSingletonClusts = sum([1 for val in clustTable.values() if val == 1])
			
		if len(OTUTable) > 0:	
			minOTUClusterSize = min(OTUTableHist.keys())
			maxOTUClusterSize = max(OTUTableHist.keys())
			aveOTUClusterSize = round(sum(OTUTable.values())/len(OTUTable),2)
			medOTUClusterSize = median(OTUTable.values())
			numOTUSingletonClusts = sum([1 for val in OTUTable.values() if val == 1])

		#if gene == "ITS":
			#pprint.pprint(clustTable)
			#input("Press Enter to continue...")  # Debugging pause
			#pprint.pprint(clustTableHist)			
			#input("Press Enter to continue...")  # Debugging pause
		numFirstOffClassifications = sum([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["Classifications"]["SH_number"] is not None and MetaDict[SampleID][gene]["Classifications"]["SH_number"].startswith("SH")])
		if gene == "RBCL":
			numFirstOffClassifications = classTable["Species"] 
		if gene == "ITS":
			SuccessfulSalvages = sum([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["Salvaging"]["Successfull"]])	
			FUNGUILD_annotations = sum([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["FUNGuild"]["Confidence_Ranking"] is not None and MetaDict[SampleID][gene]["FUNGuild"]["Confidence_Ranking"] != "NA" and MetaDict[SampleID][gene]["FUNGuild"]["Confidence_Ranking"] != "None"])
			numClusters = len(clustTable)

		print("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSummary Statistics for", gene, file=FH)
		print("\tTotal Samples:", len(MetaDict), sep="", file=FH)
		print("\tNumber of seqeuences with 2, 1, and 0 reads in the consensus:", sep="", file=FH)
		print("\t\t2 reads:", sum([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["Seqs"]["Cons"]["nSeqs"] == 2]), sep="", file=FH)
		print("\t\t1 read: ", sum([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["Seqs"]["Cons"]["nSeqs"] == 1]), sep="", file=FH)
		print("\t\t0 reads:", sum([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["Seqs"]["Cons"]["nSeqs"] not in [1, 2]]), sep="", file=FH)
		print("\tSamples with Species Level Classification (SH number):", numFirstOffClassifications, sep="", file=FH)
		if gene == "ITS" or args.cluster_plants_too:
			print("\tSamples sent through clustering:", sum([clustTable[i] for i in clustTable]), sep="", file=FH)
			print("\tSamples clustered by salvaging:", len([1 for SampleID in MetaDict if MetaDict[SampleID][gene]["Clustering"]["clustered_by_salvage"]]), sep="", file=FH)
		print("\tSamples belonging to an OTU:", sum([OTUTable[i] for i in OTUTable]), sep="", file=FH)
		print("\tSamples failed all classifications:", classTable["unclassified"], sep="", file=FH)
		if gene == "ITS":
			print("\tSamples with FUNGuild annotations:", FUNGUILD_annotations, sep="", file=FH)

		running_total = 0
		print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nClassification Metrics for", gene, file=FH)
		print("\tClassified against: ", db_used[gene], sep="", file=FH)
		running_total += classTable["Species"]
		print("\tTo Species Level:", running_total, sep="", file=FH)
		running_total += classTable["Genus"]
		print("\tTo Genus Level:", running_total, sep="", file=FH)
		running_total += classTable["Family"]
		print("\tTo Family Level:", running_total, sep="", file=FH)
		running_total += classTable["Order"]
		print("\tTo Order Level:", running_total, sep="", file=FH)
		running_total += classTable["Class"]
		print("\tTo Class Level:", running_total, sep="", file=FH)
		running_total += classTable["Phylum"]
		print("\tTo Phylum Level:", running_total, sep="", file=FH)
		running_total += classTable["Kingdom"]
		print("\tTo Kingdom Level:", running_total, sep="", file=FH)
		print("\tUnclassified:", classTable["unclassified"], sep="", file=FH)

		if gene == "ITS" or args.cluster_plants_too:
			print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nCluster Metrics for", gene, file=FH)
			print("Of the ", len(MetaDict)," total samples, ", sum([clustTable[i] for i in clustTable])," of then were clustered into ", len(clustTable), " clusters", sep="", file=FH)
			if gene =="ITS":
				print("\tand then", SuccessfulSalvages, "salvaged samples were added to those clusters:", sep=" ", file=FH)
			print("	", minClusterSize, " minimum cluster size", file=FH)
			print("	", maxClusterSize, " maximum cluster size", file=FH)
			print("	", aveClusterSize, " mean cluster size", file=FH)
			print("	", medClusterSize, " median cluster size", file=FH)
			print("	", numSingletonClusts, " number of singleton clusters", file=FH)

			print_stem_and_leaf(clustTableHist, FH)
		
		if gene == "ITS":
			print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSalvage Metrics for", gene, file=FH)
			print("Of the", len(MetaDict),"total samples grouped into", len(set(MetaDict[Sample_ID]["Metadata"]["SalvageGroup"] for Sample_ID in MetaDict)), "salvage groups:", sep=" ", file=FH)
			print("	", sum([1 for Sample_ID in MetaDict if MetaDict[Sample_ID][gene]["Salvaging"]["attempted"] == True]), " samples were suitable to attempt salvaging (i.e.: had a good enough consensus to try).", file=FH)
			print("	", SuccessfulSalvages, "samples were successfully salvaged and added to clusters.", file=FH)

		print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFinal OTU Metrics for", gene, file=FH)
		print("Of the", len(MetaDict),"total samples,", sum([OTUTable[c] for c in OTUTable]),"of then were are classified into", len(OTUTable), "OTUS:", sep=" ", file=FH)
		print("	", minOTUClusterSize, " minimum OTUs size", file=FH)
		print("	", maxOTUClusterSize, " maximum OTU size", file=FH)
		print("	", aveOTUClusterSize, " mean OTU size", file=FH)
		print("	", medOTUClusterSize, " median OTU", size, file=FH)
		print("	", numOTUSingletonClusts, " number of singleton OTUs", file=FH)

		print_stem_and_leaf(OTUTableHist, FH)

		if gene == "ITS":
			print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFUNguild Metrics for", gene, file=FH)
			print("Of the", sum([classTable[i] for i in classTable])," samples that had a taxonomy:", sep=" ", file=FH)
			print("	", funguildAttempts, " samples were sent through funguild resulting in ", fungTable["Highly Probable"] + fungTable["Probable"] + fungTable["Possible"], " successful classifications (Highly Probable + Probable + Possible).", sep="", file=FH)
			print("		", fungTable["Highly Probable"], " Highly Probable", file=FH, sep="")
			print("		", fungTable["Probable"], " Probable", file=FH, sep="")
			print("		", fungTable["Possible"], " Possible", file=FH, sep="")
			print("		", fungTable["-"], " Unclassified", file=FH, sep="")
			print("\n	The following guilds were identified:", file=FH)
			data = zip(fungTableGuild.keys(), fungTableGuild.values())
			sdata = sorted(list(data), key = lambda x: x[1], reverse=True)
			[print("		", sd[1], " ", sd[0], sep="", file=FH) for sd in sdata]
		
		print("\n\nEnd of Summary statistics for", gene, "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", file=FH)

def print_stem_and_leaf(histdict, FH):
	sl_dict = {}
	for k in sorted(histdict.keys()):
		if int(k/10) not in sl_dict:
			sl_dict[int(k/10)] = []
		[sl_dict[int(k/10)].append(k%10) for rep in range(histdict[k])]
	for k in sl_dict:
		sl_dict[k].sort()
	print("\n	Here's a stem-and-leaf plot of the distribution\n		___10s_|_1s_________________________________________________", file=FH)
	for ten in sl_dict:
		#print(ten)
		#print(sl_dict[ten])
		if len(str(ten)) == 1:#from 0 to 9
			prefix = "		 "
		elif len(str(ten)) == 2:#from 10 to 99
			prefix = "		"
		elif len(str(ten)) == 3:#from 100 to 999...
			prefix = "	   "
		elif len(str(ten)) == 4:
			prefix = "	  "
		elif len(str(ten)) == 5:
			prefix = "	 "
		
		print(prefix, ten, " | ", "".join([str(v) for v in sl_dict[ten]]), sep="", file=FH)

def write_data_summary_file(MetaDict, args, command_run, runtime):

	fileSep = "," #TAB separated file, can set this to a comma for csv.
	#IF YOU CHANGE fileSep - you must sort out FUNGuild Notes - they have commas and will mess everything up.
	end = "\n"
	if not args.onlyITS:
		end = fileSep
	OUT = os.path.join(args.output_dir, args.outputFile.rstrip("csv")+runtime+".csv")
	print("\nWriting summary file:", OUT, sep=" ", file=sys.stderr)

	version = "unknown"
	version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "VERSION")
	with open(version_file) as f:
		version = f.read().strip()

	with open(OUT, 'w') as FH:
		#first print the command so there is no confusion:
		print("#", command_run, " # version: ", version, end="\n", sep="", file=FH)
		print(fileSep.join([ #here is the column headers for the summary file
				#ITS
					#Metadata
						"Root_Code", 
						
					#consensus quality metrics
						"ITS_consensus_filename", 
						"ITS_alignQual", 
						"ITS_avQual", 
						"ITS_seqLen", 
						"ITS_nSeqs",
					
					#classification metrics - initial blast
						"ITS_blastn_filename",
						"ITS_classification_type", 
						"ITS_confidence_taxonomy", 
						"ITS_blastn_hit",
						"ITS_pident_of_top_hit", 
						"ITS_OTU",

					#clustering metrics
						"ITSx_trim_queue",
						"ITSx_trimmed",
						"ITS_clustered",
						"ITS_clustered_by_salvage",
						"ITS_cluster_ID",
						"ITS_centroid",
						"ITS_clust_size", 
						"ITS_clust_len", 
						"ITS_clust_pident", 
						"ITS_clust_strand", 
						"ITS_clust_cigar", 
						
					#Salvage metrics 
						"ITS_SalvageGroup", 
						"ITS_AttempedSalvaging", 
						"ITS_SuccessfullSalvaging",

					#FUNGuild Metrics
						"ITS_FUNGuild_Guild", 
						"ITS_FUNGuild_Confidence_Ranking" 
					]), end = end, file=FH)
		
		if not args.onlyITS: #add these in if needed
			print(fileSep.join([			
				#ITS2
					#consensus quality metrics
						"ITS2_fastq_filename", 
						"ITS2_alignQual", 
						"ITS2_avQual", 
						"ITS2_seqLen", 
						"ITS2_nSeqs",
					
					#classification metrics - initial blast
						"ITS2_blastn_filename",
						"ITS2_classification_type", 
						"ITS2_confidence_taxonomy", 
						"ITS2_blastn_hit",
						"ITS2_pident_of_top_hit",
						"ITS2_Species", 
			]),end=fileSep,file=FH)

			if args.cluster_plants_too:
				print(fileSep.join([			
					#clustering metrics
						"ITS2_clustered",#X
						"ITS2_cluster_ID",#X
						"ITS2_centroid",#X
						"ITS2_clust_size", #X
						"ITS2_clust_len", #X
						"ITS2_clust_pident", #X
						"ITS2_clust_strand", #X
						"ITS2_clust_cigar", #X	
				]),end=fileSep,file=FH)

	# #Salvage metrics 
					# 	"ITS2_SalvageGroup",#X
					# 	"ITS2_AttempedSalvaging",#X
				
			print(fileSep.join([			
				#RBCL
					#consensus quality metrics
						"RBCL_fastq_filename", 
						"RBCL_alignQual", 
						"RBCL_avQual", 
						"RBCL_seqLen", 
						"RBCL_nSeqs",
					
					#classification metrics - initial blast
						"RBCL_blastn_filename",
						"RBCL_classification_type", 
						"RBCL_confidence_taxonomy", 
						"RBCL_blastn_hit", 
						"RBCL_pident_of_top_hit",
						"RBCL_Species"
			]), end=fileSep, file=FH)

			if args.cluster_plants_too:
				print(fileSep.join([
					#clustering metrics
						"RBCL_clustered",#X
						"RBCL_cluster_ID",#X
						"RBCL_centroid",#X
						"RBCL_clust_size", #X
						"RBCL_clust_len", #X
						"RBCL_clust_pident", #X
						"RBCL_clust_strand", #X
						"RBCL_clust_cigar", #X
			]), end =fileSep, file=FH)
					# #Salvage metrics 
					# 	"RBCL_SalvageGroup", #X
					# 	"RBCL_AttempedSalvaging", #X
		if not args.onlyITS: #the final newline needs to be present if we are includeing the other ones, but it's already there for ITSonly
			print("", file=FH)		

		for SampleID in MetaDict:
			#print(MetaDict[SampleID]["ITS"]["Classifications"]["best_hit"], file=sys.stderr)
			#breaking these into pieces for ease:
			
			#general metadata: SampleID
			print(
				SampleID, sep=fileSep, end=fileSep, file=FH)
				
			#ITS consensus things: 
			print(
				MetaDict[SampleID]["ITS"]["Seqs"]["Cons"]["fasta"], 
				MetaDict[SampleID]["ITS"]["Seqs"]["Cons"]["alignmentQual"], 
				MetaDict[SampleID]["ITS"]["Seqs"]["Cons"]["avQual"], 
				MetaDict[SampleID]["ITS"]["Seqs"]["Cons"]["seqLen"],
				MetaDict[SampleID]["ITS"]["Seqs"]["Cons"]["nSeqs"], sep=fileSep, end=fileSep, file=FH)
				
			#ITS classification things: ITS_classified,			ITS_method,													ITS_SH,											ITS_taxonomy,												ITS_pident,	ITS_matchLen, ITS_bitscore,
			print(
				MetaDict[SampleID]["ITS"]["Classifications"]["blastfile"], 
				MetaDict[SampleID]["ITS"]["Classifications"]["classification_type"], 
				MetaDict[SampleID]["ITS"]["Classifications"]["taxonomy"], 
				MetaDict[SampleID]["ITS"]["Classifications"]["best_hit"],
				MetaDict[SampleID]["ITS"]["Classifications"]["best_hit"].split("\t")[2] if MetaDict[SampleID]["ITS"]["Classifications"]["best_hit"] is not None and len(MetaDict[SampleID]["ITS"]["Classifications"]["best_hit"].split("\t")) > 2 else "None", 
				MetaDict[SampleID]["ITS"]["Classifications"]["SH_number"], sep=fileSep, end=fileSep, file=FH)
				
			#ITS clustering things: RBCL_clust_ID,RBCL_clust_centroid,RBCL_clust_cigar,RBCL_clust_pident,RBCL_clust_len,RBCL_clust_size,
			print(
				MetaDict[SampleID]["ITS"]["Clustering"]["ITSx_trim_queue"],
				MetaDict[SampleID]["ITS"]["Clustering"]["ITSx_trimmed"],
				MetaDict[SampleID]["ITS"]["Clustering"]["Clustered"],
				MetaDict[SampleID]["ITS"]["Clustering"]["clustered_by_salvage"],
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["clusterID"], 
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["centroid"], 
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["clustSize"],
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["length"], 
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["pident"], 
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["strand"], 
				MetaDict[SampleID]["ITS"]["Clustering"]["cluster_data"]["cigar"], sep=fileSep, end=fileSep, file=FH)
				
			#ITS salvaged things
			print(
				MetaDict[SampleID]["Metadata"]["SalvageGroup"],
				MetaDict[SampleID]["ITS"]["Salvaging"]["attempted"],
				MetaDict[SampleID]["ITS"]["Salvaging"]["Successfull"], sep=fileSep, end=fileSep, file=FH)
			
			#ITS FUNGuild data
			for t in "Guild","Confidence_Ranking":
				print(MetaDict[SampleID]["ITS"]["FUNGuild"][t], end = fileSep, sep=fileSep, file=FH)
				
			if args.onlyITS:
				print("", end="\n", file=FH)
			else:
				#ITS2 consensus things: 
				print(
					MetaDict[SampleID]["ITS2"]["Seqs"]["Cons"]["fastq"], 
					MetaDict[SampleID]["ITS2"]["Seqs"]["Cons"]["alignmentQual"], 
					MetaDict[SampleID]["ITS2"]["Seqs"]["Cons"]["avQual"], 
					MetaDict[SampleID]["ITS2"]["Seqs"]["Cons"]["seqLen"],
					MetaDict[SampleID]["ITS2"]["Seqs"]["Cons"]["nSeqs"], sep=fileSep, end=fileSep, file=FH)
					
				#ITS2 classification things: ITS2_classified,			ITS2_method,													ITS2_SH,											ITS2_taxonomy,												ITS2_pident,	ITS2_matchLen, ITS2_bitscore,
				print(
					MetaDict[SampleID]["ITS2"]["Classifications"]["blastfile"], 
					MetaDict[SampleID]["ITS2"]["Classifications"]["classification_type"], 
					MetaDict[SampleID]["ITS2"]["Classifications"]["taxonomy"], 			
					MetaDict[SampleID]["ITS2"]["Classifications"]["best_hit"],
					MetaDict[SampleID]["ITS2"]["Classifications"]["best_hit"].split("\t")[2] if MetaDict[SampleID]["ITS2"]["Classifications"]["best_hit"] is not None and len(MetaDict[SampleID]["ITS2"]["Classifications"]["best_hit"].split("\t")) > 2 else "None", 
					MetaDict[SampleID]["ITS2"]["Classifications"]["SH_number"], sep=fileSep, end=fileSep, file=FH)
					
				# #ITS2 clustering things: RBCL_clust_ID,RBCL_clust_centroid,RBCL_clust_cigar,RBCL_clust_pident,RBCL_clust_len,RBCL_clust_size,
				if args.cluster_plants_too:
					print(
						MetaDict[SampleID]["ITS2"]["Clustering"]["Clustered"],
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["clusterID"], 
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["centroid"], 
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["clustSize"],
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["length"], 
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["pident"], 
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["strand"], 
						MetaDict[SampleID]["ITS2"]["Clustering"]["cluster_data"]["cigar"], sep=fileSep, end=fileSep, file=FH)
					
				# #ITS2 salvaged things
				# print(
				# 	MetaDict[SampleID]["Metadata"]["SalvageGroup"],
				# 	MetaDict[SampleID]["Salvaging"]["attempted"]["ITS2"], sep=fileSep, end=fileSep, file=FH)
				
				#RBCL consensus things: 
				print(
					MetaDict[SampleID]["RBCL"]["Seqs"]["Cons"]["fastq"], 
					MetaDict[SampleID]["RBCL"]["Seqs"]["Cons"]["alignmentQual"], 
					MetaDict[SampleID]["RBCL"]["Seqs"]["Cons"]["avQual"], 
					MetaDict[SampleID]["RBCL"]["Seqs"]["Cons"]["seqLen"],
					MetaDict[SampleID]["RBCL"]["Seqs"]["Cons"]["nSeqs"], sep=fileSep, end=fileSep, file=FH)
					
				#RBCL classification things: RBCL_classified,			RBCL_method,													RBCL_SH,											RBCL_taxonomy,												RBCL_pident,	RBCL_matchLen, RBCL_bitscore,
				#print(MetaDict[SampleID]["RBCL"]["Classifications"]["best_hit"], file=sys.stderr)
				print(
					MetaDict[SampleID]["RBCL"]["Classifications"]["blastfile"], 
					MetaDict[SampleID]["RBCL"]["Classifications"]["classification_type"], 
					MetaDict[SampleID]["RBCL"]["Classifications"]["taxonomy"], 
					MetaDict[SampleID]["RBCL"]["Classifications"]["best_hit"],
					MetaDict[SampleID]["RBCL"]["Classifications"]["best_hit"].split("\t")[2] if MetaDict[SampleID]["RBCL"]["Classifications"]["best_hit"] is not None and len(MetaDict[SampleID]["RBCL"]["Classifications"]["best_hit"].split("\t")) > 2 else "None", 
					MetaDict[SampleID]["RBCL"]["Classifications"]["SH_number"], sep=fileSep, end=fileSep, file=FH)
					
				#RBCL clustering things: RBCL_clust_ID,RBCL_clust_centroid,RBCL_clust_cigar,RBCL_clust_pident,RBCL_clust_len,RBCL_clust_size,
				if args.cluster_plants_too:
					print(
						MetaDict[SampleID]["RBCL"]["Clustering"]["Clustered"],
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["clusterID"], 
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["centroid"], 
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["clustSize"],
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["length"], 
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["pident"], 
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["strand"], 
						MetaDict[SampleID]["RBCL"]["Clustering"]["cluster_data"]["cigar"], sep=fileSep, end=fileSep, file=FH)
					
				# #RBCL salvaged things
				# print(
				# 	MetaDict[SampleID]["Metadata"]["SalvageGroup"],
				# 	MetaDict[SampleID]["Salvaging"]["attempted"]["RBCL"], sep=fileSep, end="\n", file=FH)
				print("", file=FH)  # End of SampleID line

def write_sequence_stats(MetaDict, args, sequence_stats_file):
	#this should read in the file at sequence_stats and then add in a column of the nSeqs for each Sample_ID and write it back out.
	with open(sequence_stats_file, "r") as seqFH:
		lines = seqFH.readlines()
	with open(sequence_stats_file, "w") as seqFH:
		#write the header line
		print(lines[0].rstrip("\n"), "nSeqs", sep="	", file=seqFH)
		for line in lines[1:]:
			Sample_ID = line.split("	")[0]
			gene = line.split("	")[1]
			if Sample_ID in MetaDict:
				nSeqs = MetaDict[Sample_ID][gene]["Seqs"]["Cons"]["nSeqs"]
				print(line.rstrip("\n"), nSeqs, sep="	", file=seqFH)
			else:
				print(line.rstrip("\n"), "NA", sep="	", file=seqFH)