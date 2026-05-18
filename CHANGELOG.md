# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [dev]

### Added

### Changed 
- updated output file to report fastq instead of fasta for ITS reads. 
- dirdepth can now parse by underscore, dash, or period. goes in order, if there's an udnerscore it uses only that, elif there's a dash, it uses only that, elif there's a period it uses only that. else it doesn't split. 

### Removed

## [2.3.0] 2026-01-21 

### Added
- added a PR checklist template
- a way to begin the pipeline at the fastq stage instead of the ab1 stage. 
- updated the readme to describe the new functionality
- added citation to the README.md
- added a CITATION.cff file
- added a flat length trimmer for raw ab1 files. 


### Changed
- there was a bug with the qual trimming, the stddev was added rather than subtracted. fixed now.
- a bug with --stddev that was forcing it to be an int if provided at the command line, not a float. fixed now. 
- updated the help page to make it clearer how setting the filters is stricter or more lenient.
- the default --prefix setting was ./ and has been changed to ""
- exposed the min length for salvaging parameter (default=150) at the commmand line
- fixed a bug with two new lines if --onlyITS is called.
- changed the Repo name to MycorrhizaTracer to match the pipeline name
- ab1s without a qual score were throwing an error, now handled appropriately. 
- simplified the conda environment to only include the primary programs and let conda sort out the dependcy tree - that should make it more portable across platforms.
- updated teh error checking for databases required, now allows less than all three databases to be specified at teh comand line. 
- was some errors with pathways depending on how a path was specified in the metadata - updated the README to address it

### Removed

## [2.2.0] 2026-01-19

### Added
- added in command line input for specifying the output files for the BOLD databases in get_BOLD_database.py
- added in a bit of test data for users to test their install with. Includes sanger reads and a metadata sheet.

### Changed
- renamed environment.yml to environment.yaml to be consistent
- the VERSION file was only accessible if the script was called from within the script dir, it is now accessible regardless of where in the computer the script is called from
- removed pathways from defaults for databases etc

### Removed

## [2.1.3] 2026-01-06

### Added

### Changed
- The name of the pipeline is now MycorrhizaTracer instead of RhizaTrace
- The name of the env is now MycorrhizaTracerEnv instead of RhizaTraceEnv

### Removed

## [2.1.2] 2025-12 First public release of pipline.

### Added

### Changed

### Removed

## [2.1.1] 2025-10-03

### Added

### Changed
- bug fix on errors that were happening for samples that only had a single sequence - those are now being handled correctly. 
- a knock-on effect of that bug fix is another bug when there are no samples to be sent through ITSx or sent through clustering. that is also fixed ans should be working now.
- updated calling pident-related flags: pidents now always specificed as percent: 0-100 (not a fraction 0-1). 

### Removed

## [2.1.0] 2025-09-02

### Added
- --min_read_len is now a command-line option to control the length of the reads that are allowed to go into the alignemnt. 
- now includes a way to output the summary stats of the raw sequencing, length and quality of each chromatogram (post filtering) and another of the nSeqs of the consensuses. 
- now outputs the counts of consenses with 2,1,and 0 reads
- now includes a test for seq records that are missing quality scores in the trimming step. (some old sanger sequencers didn't report quality scores)
- added a firstNsamples flag to run just the first N samples for troubleshooting to speed things up. This will do odd things and should never be used for final data analysis.
- added a new script that pulls the BOLD database down for the sequences described in: Jones, Laura, Alex D. Twyford, Col R. Ford, et al. 2021. ‘Barcode UK: A Complete DNA Barcoding Resource for the Flowering Plants and Conifers of the United Kingdom’. Molecular Ecology Resources 21 (6): 2050–62. https://doi.org/10.1111/1755-0998.13388.
- RhizaTraceEnv has been added in environment.yaml, it's a conda env with all the necessary programs installed: itsx, blastn, vsearch, etc, except for FUNGuild which will require a git clone call. 
- the summary output now says what classification database was used
- there is now a flag that allows one to cluster the plant reads as well as the fungi ones. 

### Changed
- changed the number of threads for the blast to 1 since there's just a single sequence ever being blasted - probably more effort to break that into args.cpus and keep track of it all than to just do a single blast. (args.cpus still used for itsx and other bits that need multiple cpus). 
- changed the blasts so that they will be submitted concurrently up to args.cpus: still all one cpu apiece, but send a couple out at a time. 
- fixed an issue with the pairwise aligner where the sequence needed to be extracted from the record. 
- the summary output now gieves better data about the salvaging.
- the salvaging is now doen against the entire set of SH numbers and centroids that were found in the samples, not just hte centroids. 
- changed the default ITS2 and RBCL database to the BOLD versions. 
- fixed a bug that failed to remove a newline in the salvage blast hits. 
- the column called ITS_SH_number is now ITS_OTU
- the command that was run now appears in the summary as well as the output
- the SequenceStats.txt file is now tab-separated instead of comma-separated.


### Removed
- removed the definitions that dealt with the odd formatting of the RBCL databases. Headers now need to match the UNITE format in all cases.

## [2.0.1] 2025-07-15 

### Added
- A file called "VERSION" which is read and written into the output file so that I know which version was used along with the command. I need to remember to update VERSION each time I push a new commit. The text in version should be strictly "2.0.0" or "2.0.0dev" etc

### Changed
- title headers were not correct for ITS2 and RBCL in the output file. That should be sorted now. Also the added in some output columns for those genes that were expected.
- species was not being pulled for RBCL as it's not the UNITE database and formated differently. That's fixed now.
- filtering for species lists is now robust to different colum seperators in the allowed species list.
- summary output for taxonomic assignments now sums each taxonomic level, so that everything classed to species is also counted when I report all the things classed to genus. previously only samples that stopped at each level were reported.

### Removed


## [2.0.0] 2025-07-14
### Notes
Many things have changed between v1.4.0 and v2.0.0. The entire code was re-written from the ground up. I pulled a couple definitions from before, but nearly everything is changed.

### Added
- Now there's an ITSx trimming step after blasting and before clustering.
- you can provide a list of speices and the blast database will be pre-emptively searhced for those hits, then the blasting will only happen against this new targeted database. Useful for incorporating regional species checklists and only looking for taxa among a restricted set. 

### Changed
- Nearly everything.
- overall logic: All samples blasted against the database, then ones that don't have a species-level hit are trimmed with ITSx and those that are successful are clustered. everything that doesn't have a species-level hit and everything that failed the ITSx trimming are blasted against the cluster centroids and joined into the best cluster they hit. The focus is on finding the best blast hit and making clusters and those two things are now the level that diversity in sites can be calculated with.
- Consensus calling: no more using Tracy, instead biopython and custom aligners and trimmers. The trimmers now trim on quality and peak height. I intend to add in one that aligns the primers and gets rid of any read-through. Trimming is based on a sliding window. Also has a check for whether the output file has been created, if so, it loads that instead of reads and trims the ab1.
- Classify: blast call is similar, but the decision making process about which read is "best" has all be re-written. Ditto for the taxonomy.
- Clustering: clusters with vsearch as before, but only uses the best of the best sequences to build clusters from, then adds unclustered sequences in via blast to the centroids during salvaging.
- Salvaging: salvaging is now done for all samples that weren't clustered. If a sample finds a cluster and the sample had no good blast hit, the classification of the sample is pulled from the centroid, but if the sample has a good hit then the sample's original hit is used. 
- FUNGuild: mostly similar, now outputs only a couple of the fields though
- Output file: lots of fields have been removed in an effor to keep the output tidy. Some have changed form i.e.: no more "primary"/"secondary" etc, but now it's the number of sequences in the consensus: 2, 1, or 0.
- sample input sheet: got rid of a bunch of columns that weren't necessary.

### Removed
-sample checking before the pipeline begins - currently left to the user to check that the samples are in order. I want to add this back in eventually though.

## [1.4.0] - 2025-07-05
### Notes
- this is the final release of version 1. 

### Added
- added this change log
- log files now have the datetime the run started

### Changed
-default tracy trimming stringency of BestReads changed from 2 to 3. 

### Removed


## [1.3.5] - 2025-01-13

### Added

### Changed

### Removed
