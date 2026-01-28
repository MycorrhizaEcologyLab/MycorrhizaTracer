#!/usr/bin/env python3
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description="Extract the BOLD database sequences for UK plants for ITS2 and RBCL and write them to a FASTA file each.")
	#general arguments
	parser.add_argument("-i", "--ITS2_output", type=str, required=False, default = "./BOLD_database_ITS2.fasta", help="Path to the ITS2 output file. Default is %(default)s.")
	parser.add_argument("-r", "--RBCL_output", type=str, required=False, default = "./BOLD_database_RBCL.fasta", help="Path to the RBCL output file. Default is %(default)s.")
	return parser.parse_args()

args = parse_args()

#This is a bit of code to query the BOLD API for sequences of specific samples and write them to FASTA files.
#it's hard-coded for the specific samples that are described in the following paper which are UK species of plants:
#   Jones, Laura, Alex D. Twyford, Col R. Ford, et al. 2021. 
#   ‘Barcode UK: A Complete DNA Barcoding Resource for the Flowering Plants and Conifers of the United Kingdom’. 
#   Molecular Ecology Resources 21 (6): 2050–62. https://doi.org/10.1111/1755-0998.13388.


#need to use python to query the BOLD api for FPUK001-14 to FPUK1362-20; POWNA001-10to POWNA3220-13; POWNB001-10 to POWNB237-10. 
#http://v3.boldsystems.org/index.php/API_Public/sequence?ids=POWNB001-10|

import requests
import os
import sys
import time
import xml.etree.ElementTree as ET

#open a file called BOLD_database_RBCL.fasta and BOLD_database_ITS2.fasta in the current directory
#query the BOLD API for the sequences of the samples in the file
#check each download and write to the database if RBCL is in the sequence header or to the other databse if ITS2 is in the sequence header
def get_BOLD_database():
	# Define the BOLD API URL
	api_url = "http://v3.boldsystems.org/index.php/API_Public/combined" 
	# Define the output files
	rbcl_output_file = args.RBCL_output
	its2_output_file = args.ITS2_output
	# Define the sample ranges and years
	sample_ranges = {"FPUK": 1362, "POWNA": 3220, "POWNB": 237}

	sample_years = {"FPUK": range(14, 21), "POWNA": range(10, 14), "POWNB": range(10, 11)}  # Years for each sample type

	with open(rbcl_output_file, 'w') as rbcl_file, open(its2_output_file, 'w') as its2_file:
		BATCH_SIZE = 50  # You can adjust this
		for name, max_id in sample_ranges.items():
			for year in sample_years[name]:
				# Use the known max only for the final year, otherwise use 9999
				if year == max(sample_years[name]):
					upper = max_id
				else:
					upper = 9999
				no_data_count = 0
				sample_ids = [f"{name}{i:03d}-{year}" for i in range(1, upper + 1)]
				for batch_start in range(0, len(sample_ids), BATCH_SIZE):
					batch = sample_ids[batch_start:batch_start+BATCH_SIZE]
					ids_str = "|".join(batch)
					api_query = f"{api_url}?ids={ids_str}"
					print(f"Processing batch: {batch[0]} ... {batch[-1]}", file=sys.stderr)
					# Make the API request
					response = requests.get(api_query)
					time.sleep(0.5)
					# Check if the request was successful
					if response.status_code == 200:
						# Parse the response as a multiline-fasta:
						# Assuming the response is in FASTA format, split by '>' to get individual sequences	
						if response.text.strip() == "":
							print(f"\tNo data found for batch {batch[0]} ... {batch[-1]}. Skipping...")
							no_data_count += 1
							if no_data_count >= 50:
								print(f"\t50 consecutive no data for {name}-{year}. Assuming end of samples for this year.", file=sys.stderr)
								break
							continue
						no_data_count = 0  # Reset counter if data is found
						if response.text.strip().startswith("<?xml"):
							fasta_records = bold_xml_to_unite_fasta(response.text)
							#print(fasta_records, file=sys.stderr)
							for rec in fasta_records:
								# You can filter for RBCL/ITS2 by marker code if needed
								if "rbcl" in rec.lower():
									rbcl_file.write(rec + "\n")
								elif "its2" in rec.lower():
									its2_file.write(rec + "\n")
						# else:#this really only happens if the API returns a FASTA format which only happens when the query is a sequence api query not a combined api query, this will probably not be run. 
						#	 # Split the response into individual sequences
						#	 fasta_sequences = response.text.strip().split('>')[1:]  # Skip the first empty split
						#	 for fasta_sequence in fasta_sequences:
						#		 # Split the sequence into header and data
						#		 fasta_parts = fasta_sequence.split('\n', 1)
						#		 if len(fasta_parts) < 2:
						#			 print(f"\tInvalid FASTA format for one of the entries. Skipping...")
						#			 continue
						#		 fasta_header = fasta_parts[0].strip()
						#		 fasta_data = fasta_parts[1].strip()
						#		 # Check if the header contains RBCL or ITS2
						#		 if "rbcl" in fasta_header.lower():  
						#			 rbcl_file.write(f">{fasta_header}\n{fasta_data}\n")
						#		 elif "its2" in fasta_header.lower():
						#			 its2_file.write(f">{fasta_header}\n{fasta_data}\n")
					else:
						print(f"\tHTTP error {response.status_code} for batch {batch[0]} ... {batch[-1]}", file=sys.stderr)
						print(f"\tapi call: {api_query}", file=sys.stderr)


def bold_xml_to_unite_fasta(xml_text):
	records = []
	root = ET.fromstring(xml_text)
	for record in root.findall('record'):
		# Extract taxonomy
		try:
			genus = record.find('./taxonomy/genus/taxon/name').text
			species = record.find('./taxonomy/species/taxon/name').text
		except AttributeError:
			continue  # Skip if genus or species missing

		# Format genus_species
		genus_species = f"{species.replace(' ', '_')}"
		# Extract process/sample ID
		sample_id = record.find('processid').text if record.find('processid') is not None else "NA"
		# Extract taxonomy fields
		kingdom = record.find('./taxonomy/kingdom/taxon/name')
		phylum = record.find('./taxonomy/phylum/taxon/name')
		class_ = record.find('./taxonomy/class/taxon/name')
		order = record.find('./taxonomy/order/taxon/name')
		family = record.find('./taxonomy/family/taxon/name')
		# Build taxonomy string
		tax_str = f"k__{kingdom.text if kingdom is not None else 'Viridiplantae'};" \
				  f"p__{phylum.text if phylum is not None else 'NA'};" \
				  f"c__{class_.text if class_ is not None else 'NA'};" \
				  f"o__{order.text if order is not None else 'NA'};" \
				  f"f__{family.text if family is not None else 'NA'};" \
				  f"g__{genus};s__{species.replace(' ', '_')}"
		# Extract sequence(s)
		for seq in record.findall('./sequences/sequence'):
			nucleotides = seq.find('nucleotides')
			markercode = seq.find('markercode')
			gene = markercode.text if markercode is not None else "NA"
			if nucleotides is not None and nucleotides.text:
				header = f">{genus_species}|NA|{sample_id}|{gene}|{tax_str}"
				sequence = nucleotides.text.replace('\n', '').replace("-", "")
				records.append(f"{header}\n{sequence}")
	return records



if __name__ == "__main__":
	get_BOLD_database()