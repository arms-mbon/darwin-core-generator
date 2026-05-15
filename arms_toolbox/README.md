## Description 

Here you find the code in arms_toolbox modified to run specifically on the ARMS-MBON data, 18S, ITS, and COI as provided by the PEMA pipeline, more specifically PEMA v 2.1.4.  This means it is turned to the column names and file formatting (OTU/ASV, taxonomy tables, and fasta files) from this version of PEMA, and is tuned to the names we have given to the respective files. These files can be found in  [analysis_release_001](https://github.com/arms-mbon/analysis_release_001) and [data_release_001](https://github.com/arms-mbon/data_release_001) for the first data release, and [analysis_release_002](https://github.com/arms-mbon/analysis_release_002) and [data_release_002](https://github.com/arms-mbon/data_release_002) for the second. The code takes the data from these repositories, and additionally some data is hard coded in the input file main.py and in the . 

We have applied some checks to the taxonomic names provided by PEMA, to confirm to the expectations of EurOBIS.
* We only take results that have a quantity value of >1
* For COI, where confidence levels are reported for each name in the taxonomic classification assigned to the ASVs, we only report the ScientificName in the occurrence.csv file that has a level >=0.8. So if the ScientificName found by PEMA is less than this, we go down the taxonomic classification (i.e. to the next less-specific level) and check that one. The first that has a level >=0.8 is the one accepted. 
* Where a scientific name is not also found in the Catalogue of Life, then we (again, for COI) go down the taxonomic classification (i.e. to the next less-specific level) and check that one. The first name that _is_ in the CoF is accepted as the ScientificName to report
* Where a ScientificName cannot be found in WoRMS (there is no aphia ID), then we simply return no scientificNameID
* In the occurrence.csv file, the full taxonomic classiciation returned by PEMA is included, as well as the NCBI ID for the scientific name that it found. So this original information is preserved. 
* In the emof.csv file, for COI we report the original scientific name and confidence level that was found by PEMA. So this original information is preserved.
We use APIs from WoRMS, CoL, and ENA to get the various IDs that are required when doing this if-then-elseing. 

Other notes
* The emof schema can be optimised to the gene type.
* This code is designed to be run over several input files whoes name differ by gene type (18S and ITS) and date (e.g. July2020).
* The output is added to data/outputs.
* Checks for duplicates and missing information are performed - for duplicates, a warning is logged and the first only is taken.
* We only take those results with a "quantity" value in the OTU table that is >1.
* A check on ENA to find the sample accession numbers for the run accession numbers that we have in the omics input is done.
* ITS, 18S, and COI differ in that (1) the input files come from different repo folders, (2) COI has more in the emof and takes that from a PEMA file that only exists for this one

Please read the comments in the code before using it!

## The IDs in the input files used for this data_release_001 processing:

Extended_final_table (col 1)
* For 18S and ITS, it is Otu###
* For COI it is ASV_##:#########(long number)

Fasta files (these are read in via main.py)
* For 18S and ITS it is the same Otu###
* for COI it is (grrrrr!) #######(long number)_####, i.e. the second part of that from the Extended_final_table with a _(number) appended

Tax_assigments
* For 18S and ITS there is no such file
* For COI it is the same as in the fasta file

## Corrections applied in the COI dataset of data_release_001

Since the code was not working properly, there were species level identifications in the dataset with very low confidence estimates. 
These identifications have now been corrected using [ARMS_issues_COI_2018_2020.R](https://github.com/arms-mbon/darwin-core-generator/blob/main/arms_toolbox/ARMS_issues_COI_2018_2020.R). The image is also available [here](https://github.com/arms-mbon/darwin-core-generator/blob/main/arms_toolbox/ARMS_issues_COI_2018_2020.RData). 

Also, in cases where confidence estimates were higher than 0.8 only at the Kingdom level (i.e. Eukaryota), these occurrences are now reported as Biota incertae sedis because it may be very likely that they actually belong in another Kingdom and not even in the Eukaryota kingdom.

The updated dataset can be found [here](https://github.com/arms-mbon/data_release_001/blob/main/ARMS_COI_Occurrence_corrected.txt).

After the confidence estimate corrections, the dataset was uploaded in the [geographic outlier detection tool of OBIS](https://ednaqc.obis.org).
Based on the tool, density values close to zero indicate occurrence of taxa outside of their known distributions. 
Also, suitability values close to zero indicate occurrence of taxa outside of their known temperature preference. 
Low density scores combined with high suitability scores may indicate a detection of an introduced taxon. However, low density scores combined with low suitability scores indicate likely erroenous identifications. Based on the above, annotations were either accepted or rejected and relevant information was added in the comments. The annotated are provided [here](https://github.com/arms-mbon/data_release_001/blob/main/ARMS_COI_Occurrence_corrected_annotations.json) and [here](https://github.com/arms-mbon/data_release_001/blob/main/ARMS_COI_Occurrence_corrected_annotations_old_format.json).

## Corrections applied in the COI dataset of data_release_002

As mentioned in the previous section, the COI dataset for data_release_002 has been corrected, using [ARMS_issues_COI_2020_2021.R](https://github.com/arms-mbon/darwin-core-generator/blob/main/arms_toolbox/ARMS_issues_COI_2020_2021.R). The image is also available [here](https://github.com/arms-mbon/darwin-core-generator/blob/main/arms_toolbox/ARMS_issues_COI_2020_2021.RData).

The updated dataset can be found [here](https://github.com/arms-mbon/data_release_002/blob/main/ARMS_COI_Occurrence_corrected.txt).


