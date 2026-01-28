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

#../Classify_Sanger_OTUs.py --metadata ./Sanger_database_demoset.csv --prefix ../ --cpus 6 --ITS_db ./dbs/sh_general_release_dynamic_s_04.04.2024.fasta --ITS2_db ./dbs/BOLD_database_ITS2.fasta --RBCL_db ./dbs/BOLD_database_RBCL.fasta --salvage --run_funguild # version: 2.2.0


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary Statistics for ITS
	Total Samples:31
	Number of seqeuences with 2, 1, and 0 reads in the consensus:
		2 reads:20
		1 read: 11
		0 reads:0
	Samples with Species Level Classification (SH number):24
	Samples sent through clustering:6
	Samples clustered by salvaging:5
	Samples belonging to an OTU:31
	Samples failed all classifications:0
	Samples with FUNGuild annotations:23

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classification Metrics for ITS
	Classified against: sh_general_release_dynamic_s_04.04.2024.fasta
	To Species Level:24
	To Genus Level:25
	To Family Level:28
	To Order Level:31
	To Class Level:31
	To Phylum Level:31
	To Kingdom Level:31
	Unclassified:0

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cluster Metrics for ITS
Of the 31 total samples, 6 of then were clustered into 2 clusters
	and then 5 salvaged samples were added to those clusters:
	 1  minimum cluster size
	 5  maximum cluster size
	 3.0  mean cluster size
	 3.0  median cluster size
	 1  number of singleton clusters

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 15

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Salvage Metrics for ITS
Of the 31 total samples grouped into 11 salvage groups:
	 5  samples were suitable to attempt salvaging (i.e.: had a good enough consensus to try).
	 5 samples were successfully salvaged and added to clusters.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Final OTU Metrics for ITS
Of the 31 total samples, 31 of then were are classified into 7 OTUS:
	 1  minimum OTUs size
	 6  maximum OTU size
	 4.43  mean OTU size
	 5  median OTU 4
	 1  number of singleton OTUs

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 1445566

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNguild Metrics for ITS
Of the 31  samples that had a taxonomy:
	23 samples were sent through funguild resulting in 23 successful classifications (Highly Probable + Probable + Possible).
		15 Highly Probable
		8 Probable
		0 Possible
		0 Unclassified

	The following guilds were identified:
		15 |Ectomycorrhizal|
		5 Plant PathogenNone|Plant Saprotroph|NoneUndefined SaprotrophNoneWood Saprotroph
		2 Plant SaprotrophNoneWood Saprotroph
		1 EctomycorrhizalNone|Endophyte|NoneEricoid MycorrhizalNoneOrchid MycorrhizalNonePlant SaprotrophNoneUndefined Saprotroph


End of Summary statistics for ITS 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary Statistics for ITS2
	Total Samples:31
	Number of seqeuences with 2, 1, and 0 reads in the consensus:
		2 reads:23
		1 read: 4
		0 reads:4
	Samples with Species Level Classification (SH number):0
	Samples belonging to an OTU:14
	Samples failed all classifications:9

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classification Metrics for ITS2
	Classified against: BOLD_database_ITS2.fasta
	To Species Level:14
	To Genus Level:16
	To Family Level:17
	To Order Level:18
	To Class Level:22
	To Phylum Level:22
	To Kingdom Level:22
	Unclassified:9

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Final OTU Metrics for ITS2
Of the 31 total samples, 14 of then were are classified into 2 OTUS:
	 3  minimum OTUs size
	 11  maximum OTU size
	 7.0  mean OTU size
	 7.0  median OTU 3
	 0  number of singleton OTUs

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 3
		 1 | 1


End of Summary statistics for ITS2 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary Statistics for RBCL
	Total Samples:31
	Number of seqeuences with 2, 1, and 0 reads in the consensus:
		2 reads:21
		1 read: 5
		0 reads:5
	Samples with Species Level Classification (SH number):12
	Samples belonging to an OTU:12
	Samples failed all classifications:7

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classification Metrics for RBCL
	Classified against: BOLD_database_RBCL.fasta
	To Species Level:12
	To Genus Level:21
	To Family Level:23
	To Order Level:24
	To Class Level:24
	To Phylum Level:24
	To Kingdom Level:24
	Unclassified:7

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Final OTU Metrics for RBCL
Of the 31 total samples, 12 of then were are classified into 2 OTUS:
	 3  minimum OTUs size
	 9  maximum OTU size
	 6.0  mean OTU size
	 6.0  median OTU 3
	 0  number of singleton OTUs

	Here's a stem-and-leaf plot of the distribution
		___10s_|_1s_________________________________________________
		 0 | 39


End of Summary statistics for RBCL 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

```
