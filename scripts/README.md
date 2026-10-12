## Description 

### The xxx.RData files
Are the R scripts (see below) in R format

 ### The xxx.R files
This code is _not_ part of the DwC generator, however it _is_ part of the ARMS-MBON struggle to create DwC files for submission to EurOBIS with suitable taxonomic information. These scripts were created to make modifications to the taxonomic names arising out of the ARMS-MBON COI data as provided by the PEMA pipeline v 2.1.4. These corrections were at the same time also implemented in the DwC generator code run for data release 002; further modifications - hopefully the final ones - will be in place for data release 003, and these will be retrospectively applied to data releases 001 and 002. Hence, the code here is stand-alone and was not as-is included in the DwC generator code. The file and column names are as produced by ARMS-MBON using that version of PEMA, and there are some hard-coded taxonomic name corrections included. 


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

The COI dataset for data_release_002 has been corrected, using [ARMS_issues_COI_2020_2021.R](https://github.com/arms-mbon/darwin-core-generator/blob/main/arms_toolbox/ARMS_issues_COI_2020_2021.R). The image is also available [here](https://github.com/arms-mbon/darwin-core-generator/blob/main/arms_toolbox/ARMS_issues_COI_2020_2021.RData).

The updated dataset can be found [here](https://github.com/arms-mbon/data_release_002/blob/main/ARMS_COI_Occurrence_corrected.txt).

## Corrections applied in the 18S dataset of data_release_001

The dataset was corrected since the primer references were wrong. 

The updated dataset can be found [here](https://github.com/arms-mbon/data_release_001/blob/main/ARMS_18S_DNAextension_corrected.txt).

