## Description 

Here you find the code in arms_toolbox modified to run specifically on the ARMS data of data_release_001 - this was necessary mainly for the COI outputs of PEMA, but necessitated a few tweaks also to the 18S and ITS codes. 

Note that this code is designed to run on the outputs of PEMA v 2.1.4 and the column names and file formatting (OTU/ASV, taxonomy tables, and fasta files) therefrom, as well as the file names. Input also comes from fasta files archived in the MDA and from ARMS logsheets from which event and sampling information are extracted.

The emof schema can be optimised to the gene type.

This code is designed to be run over several input files whoes name differ by gene type (18S and ITS) and date (e.g. July2020).
The output is added to data/outputs.

Checks for duplicates and missing information are performed - for duplicates, a warning is logged and the first only is taken.

We only take those results with a "quantity" value in the OTU table that is >1.  

A check on ENA to find the sample accession numbers for the run accession numbers that we have in the omics input is done.
ITS, 18S, and COI differ in that 
* the input files come from different repo folders
* COI has more in the emof and takes that from a PEMA file that only exists for this one

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

