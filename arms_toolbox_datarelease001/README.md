Here you find the code in arms_toolbox modified to run specifically on the ARMS data of data_release_001 - this was necessary mainly for the COI outputs of PEMA, but necessitated a few tweaks also to the 18S and ITS codes. 
Note that this code is designed to run on the outputs of PEMA v 2.1.4 and the column names and file formatting (OTU/ASV, taxonomy tables, and fasta files) therefrom.

Please read the comments in the code -- there is some hard-coding of file names and file location in these codes, and we only take those results with a "quantity" value in the OTU table that is >1.  

