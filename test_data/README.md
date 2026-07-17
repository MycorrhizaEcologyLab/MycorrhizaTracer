A small test dataset for the pipeline that includes some high and some low quality sequences. 
You'll need to unzip the UNITE database in the dbs directory (it was too big for github unzipped)

Suggested test run from within the ./test_data/ directory:
```
../Classify_Sanger_OTUs.py \
    --metadata ./Sanger_database_demoset.csv \
    --prefix ../ \
    --output_dir ./TEST_OUTPUT/ \
    --cpus 6 \
    --ITS_db ./dbs/sh_general_release_dynamic_s_04.04.2024.fasta \
    --ITS2_db ./dbs/BOLD_database_ITS2.fasta \
    --RBCL_db ./dbs/BOLD_database_RBCL.fasta \
    --salvage \
    --run_funguild \
```


It should create a directory called ./TEST_OUTPUT/ with a variety of files in it. One of which is called Sanger_Summary.YYYY-MM-DD_HHMM.txt and should contain the following:



```
#../Classify_Sanger_OTUs.py --metadata ./Sanger_database_demoset.csv --prefix ../ --output_dir ./TEST_OUTPUT/ --cpus 6 --ITS_db ./dbs/sh_general_release_dynamic_s_04.04.2024.fasta --ITS2_db ./dbs/BOLD_database_ITS2.fasta --RBCL_db ./dbs/BOLD_database_RBCL.fasta --salvage --run_funguild # version: 2.3.1


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary Statistics for ITS
	Total Samples:31
	Number of seqeuences with 2, 1, and 0 reads in the consensus:
		2 reads:22
		1 read: 9
		0 reads:0
	Samples with Species Level Classification (SH number):21
	Samples sent through clustering:7
	Samples clustered by salvaging:1
	Samples belonging to an OTU:29
	Samples failed all classifications:0
	Samples with FUNGuild annotations:21

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classification Metrics for ITS
	Classified against: sh_general_release_dynamic_s_04.04.2024.fasta
	To Species Level:21
	To Genus Level:22
	To Family Level:26
	To Order Level:28
	To Class Level:30
	To Phylum Level:31
	To Kingdom Level:31
	Unclassified:0

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cluster Metrics for ITS
Of the 31 total samples, 7 of then were clustered into 4 clusters
	and then 1 salvaged samples were added to those clusters:
	 1  minimum cluster size
	 4  maximum cluster size
	 1.75  mean cluster size
	 1.0  median cluster size
	 3  number of singleton clusters

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 1114

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Salvage Metrics for ITS
Of the 31 total samples grouped into 11 salvage groups:
	 3  samples were suitable to attempt salvaging (i.e.: had a good enough consensus to try).
	 1 samples were successfully salvaged and added to clusters.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Final OTU Metrics for ITS
Of the 31 total samples, 29 of then were are classified into 9 OTUS:
	 1  minimum OTUs size
	 7  maximum OTU size
	 3.22  mean OTU size
	 3  median OTU 3
	 3  number of singleton OTUs

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 111334457

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNguild Metrics for ITS
Of the 31  samples that had a taxonomy:
	21 samples were sent through funguild resulting in 21 successful classifications (Highly Probable + Probable + Possible).
		14 Highly Probable
		7 Probable
		0 Possible
		0 Unclassified

	The following guilds were identified:
		14 |Ectomycorrhizal|
		4 Plant PathogenNone|Plant Saprotroph|NoneUndefined SaprotrophNoneWood Saprotroph
		2 Plant SaprotrophNoneWood Saprotroph
		1 EctomycorrhizalNone|Endophyte|NoneEricoid MycorrhizalNoneOrchid MycorrhizalNonePlant SaprotrophNoneUndefined Saprotroph


End of Summary statistics for ITS 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary Statistics for ITS2
	Total Samples:31
	Number of seqeuences with 2, 1, and 0 reads in the consensus:
		2 reads:22
		1 read: 5
		0 reads:4
	Samples with Species Level Classification (SH number):0
	Samples belonging to an OTU:15
	Samples failed all classifications:9

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classification Metrics for ITS2
	Classified against: BOLD_database_ITS2.fasta
	To Species Level:15
	To Genus Level:17
	To Family Level:18
	To Order Level:20
	To Class Level:22
	To Phylum Level:22
	To Kingdom Level:22
	Unclassified:9

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Final OTU Metrics for ITS2
Of the 31 total samples, 15 of then were are classified into 3 OTUS:
	 1  minimum OTUs size
	 12  maximum OTU size
	 5.0  mean OTU size
	 2  median OTU 1
	 1  number of singleton OTUs

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 12
		 1 | 2


End of Summary statistics for ITS2 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary Statistics for RBCL
	Total Samples:31
	Number of seqeuences with 2, 1, and 0 reads in the consensus:
		2 reads:22
		1 read: 5
		0 reads:4
	Samples with Species Level Classification (SH number):14
	Samples belonging to an OTU:14
	Samples failed all classifications:6

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classification Metrics for RBCL
	Classified against: BOLD_database_RBCL.fasta
	To Species Level:14
	To Genus Level:22
	To Family Level:24
	To Order Level:25
	To Class Level:25
	To Phylum Level:25
	To Kingdom Level:25
	Unclassified:6

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Final OTU Metrics for RBCL
Of the 31 total samples, 14 of then were are classified into 2 OTUS:
	 4  minimum OTUs size
	 10  maximum OTU size
	 7.0  mean OTU size
	 7.0  median OTU 4
	 0  number of singleton OTUs

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 4
		 1 | 0


End of Summary statistics for RBCL 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


```
