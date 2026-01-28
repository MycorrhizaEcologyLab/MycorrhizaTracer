# Classify_Sanger_OTUs

This Python pipeline processes raw Sanger sequencing data (ab1 files) or manually curated fasta/fastq files, generates consensus sequences, classifies them taxonomically using BLASTN against reference databases, clusters sequences into OTUs, and optionally salvages low-quality reads by matching them to cluster centroids. The pipeline is designed for fungal (full ITS), and plant (ITS2, and RBCL) genes, but could be adapted for other loci. Functional annotation of fungi is performed using FUNGuild.

## Overview

- **Consensus calling**: Generates consensus sequences from ab1 files or uses manually curated fasta/fastq consensus files.
- **Classification**: Uses BLASTN to assign taxonomy to each consensus sequence, with configurable thresholds for species/genus/family/etc.
- **Clustering**: Clusters classified sequences into OTUs using vsearch, assigning SH numbers and consensus taxonomy.
- **Salvaging**: Optionally attempts to classify low-quality or unclassified reads by blasting them against cluster centroids.
- **Functional annotation**: Runs FUNGuild to assign ecological guilds to classified fungal sequences.
- **Comprehensive output**: Produces summary tables, per-sample files, and logs.

## Required Files

- `Classify_Sanger_OTUs.py`
- `Classify_Sanger_OTUs_definitions_module.py`
- A metadata CSV file linking unique sample IDs to ab1 file locations and sample information (see template below).

## Required Installs

- **Conda**: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
- **Python 3.11**: https://www.python.org/downloads/
- **Biopython 1.83**: https://biopython.org/wiki/Getting_Started
- **BLAST+ 2.14.0**: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata
- **vsearch 2.21.11**: https://github.com/torognes/vsearch
- **FUNGuild 1.1.0**: https://github.com/UMNFuN/FUNGuild
- **ITSx 1.1.3**: https://microbiology.se/software/itsx/
- **Reference databases**: UNITE, UNITE-all, an rbcL database, or other FASTA-formatted databases with headers as found in the UNITE general release: `Genus_species|{unused_code}|{SH_number_or_similar}|{unused_code}|k__Fungi;p__Ascomycota;...;g__Genus;s__species`.

The appropriate versions of Python, Biopython, BLAST+, vsearch, and ITSx are included in the provided conda environment file (`environment.yaml`).  
To set up the environment (only needed once), run:

```bash
conda env create -f environment.yaml
```

To access the environment prior to running the pipeline, run:
```bash
conda activate MycorrhizaTracerEnv
```

Install FUNGuild by following the instructions here: https://github.com/UMNFuN/FUNGuild



## Metadata File

The metadata CSV should include at least the following columns:

| Column Name     | Required? | Notes                                                                                  |
|-----------------|-----------|---------------------------------------------------------------------------------------|
| Sample_ID       | YES       | Unique sample identifier                                                              |
| Plot.           | optional  | Used for salvaging: restricts salvaged assignments to groups already present in plot  |
| directory       | YES       | Directory containing the ab1 files                                                    |
| ITS_For_file    | depends   | Relative path to ITS forward ab1 file                                                 |
| ITS_Rev_file    | depends   | Relative path to ITS reverse ab1 file                                                 |
| ITS2_For_file   | depends   | Relative path to ITS2 forward ab1 file                                                |
| ITS2_Rev_file   | depends   | Relative path to ITS2 reverse ab1 file                                                |
| RBCL_For_file   | depends   | Relative path to RBCL forward ab1 file                                                |
| RBCL_Rev_file   | depends   | Relative path to RBCL reverse ab1 file                                                |
| ...             | optional  | Any other metadata columns you wish to include                                        |

See the provided template "Metadata_Template.xlsx" for details, this can be filled out and will need to be saved in CSV format.
If you have manually curated fasta/fastq files: the pipeline can begin with them if a second metadata sheet in the same style is filled out. The Sample_IDs in it must appear in both metadata sheets and the plot column as well as all Rev columns should be empty. Only use For columns for manually curated consensus sequences. Pass this additional file to the pipeline with --manual_consensus_file and all ab1 quality trimming, alignment, and consensus building steps will be skipped for that sample in favour of using the curated file you provided.

## Required Flags

- `--meta` (metadata CSV file)
- `--output_dir` (directory for all outputs)
- `--prefix` (prefix for sequencing file paths)
- `--outputFile` (summary output file)
- `--ITS_db`, `--ITS2_db`, `--RBCL_db` (paths to reference databases)

## Recommended Flags

- `--resume` Re-run pipeline without overwritting previous ouputs. Delete any files you want repeated, but make sure to deleter everything downstream of them as well.
- `--adjustPidents` (force global alignments for BLASTN)
- `--dirDepth` (organize output files into subdirectories for large datasets)

## Metadata Checking

The script does not perform checks to ensure all files are accessible and metadata is consistent. Please ensure that your metadata file is free of errors. Though this may be implemented eventually.

## Consensus Calling

Consensus sequences are generated using in-built QC trimming and aligning algorithms implemented with BioSeq. If consensus quality is low or ab1 files are missing, the best available read is used. Samples with no usable sequence are discarded.

Relevant flags: `--qual`, `--window_size`, `--stddev_cutoff`, `--manual_consensus_file`

## Classification Algorithm

BLASTN is used to classify each consensus. Hits are filtered by length, coverage, and percent identity. The best hit is chosen based on percent identity and length, and consensus taxonomy is assigned at the highest supported level.

Relevant flags: `--minMatchLength`, `--differential`, `--minSCOV`, `--minQCOV`, `--min_pident_species`, `--min_pident_genus`, `--min_pident_family`, `--min_pident_phylum`, `--min_consensus_blast`, `--no_Incertae_Sedis`, `--adjustPidents`, `--species_list_ITS`, `--species_list_ITS2`, `--species_list_RBCL`

## Clustering

Sequences are clustered into OTUs using vsearch. Each cluster is assigned a consensus taxonomy and a unique SH number.

Relevant flags: `--minSeqLengthForCluster`, `--minAverageQualityForCluster`, `--percentIdentityForCluster`

## Salvaging

Unclassified or low-quality samples can be optionally salvaged by blasting them against cluster centroids. Assignments are restricted to SalvageGroups if specified.

Relevant flag: `--salvage` `--min_read_len_salvage`

## FUNGuild

Classified fungal sequences are annotated with FUNGuild. Set the path to the executable with `--FUNGuild_executable` if needed.

Relevant flags: `--FUNGuild_executable`, `--run_FUNGUild`

## Output Files

- Consensus FASTQ/FASTA files
- BLASTN output files
- Cluster files (*.uc)
- FUNGuild annotation files
- Summary tables and logs

The main output file is named with the -f or --outputFile flag (default: "Sanger_Output.csv"). The file "Output_key.README" is a description of each column that will be found therein.

## Usage

```sh
python Classify_Sanger_OTUs.py --metadata samples.csv --output_dir ./results/ --outputFile output_file.csv --ITS_db UNITE_location.fasta --ITS2_db UNITE_ALL_location.fasta --RBCL RBCL_location.fasta --salvage --run_funguild 
```

For full options, run:
```sh
python Classify_Sanger_OTUs.py -h
```

## Test Data

There is a provided dataset and metadata sheet in ./test_data/ that can be used to test the install with. See the README.txt there for a suggested call.


## Notes

- Taxonomy headers in reference FASTA files must be formatted as UNITE general fasta release
- Salvaging should always be interpreted with caution.
- For best results, review the summary and logs for warnings or errors.
- If you are adapting the pipeline for taxa that are not fungi or plants or using gene regions other than fungal ITS, you MUST use the columns in the metadata sheet assigned to ITS2 or RBCL (and also specify an appropriate --ITS2_db or --RBCL_db). This is because the chromatograms in the ITS columns will be passed to ITSx which will search for the ITS region. If it is missing or incomplete then that sample will fail. This step is skipped for the ITS2 and RBCL chromatograms, and so is appropriate for all other regions and taxa. 

---

For questions or feature requests, please contact the script author or open an issue on github


## Citation

Brekke, T.D., Weeks, T.L., Barber, R.A., Thomson, I., Gooda, R., Gargiulo, R., Delhaye, G., Andrew, C., Kowal, J., Bidartondo, M., Suz, L. M. (in prep) MycorrhizaTracer: A BIOINFORMATIC PIPELINE FOR FUNGI AND PLANT CLASSIFICATION OF SANGER DNA SEQUENCES