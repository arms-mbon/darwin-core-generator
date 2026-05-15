
# ============================================================
'R code for correction of issues in ARMS-MBON datasets

Christina Pavloudi
christina.pavloudi@embrc.eu
https://cpavloud.github.io/mysite/

	Copyright (C) 2026 Christina Pavloudi
  
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'

# =============================================================

############################LOAD LIBRARIES #######################################

# List of packages needed
.packages = c("naniar", "dplyr", "obistools", "tidyverse")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

packageVersion("naniar")
packageVersion("dplyr")
packageVersion("obistools")
packageVersion("tidyverse")

################################################################################
################################################################################
################################################################################

#getwd()
#setwd("")

################################################################################
################################################################################
################################################################################

#read input DwC files
#ARMS_COI_2020_2021_occurrence <- read.csv("dwca-arms_coi_2020-22-v1/occurrence.txt", sep = "\t", header=TRUE) 
#ARMS_COI_2020_2021_eMOF <- read.csv("dwca-arms_coi_2020-22-v1/extendedmeasurementorfact.txt", sep = "\t", header=TRUE) 
ARMS_COI_2018_2020_occurrence <- read.csv("dwca-arms_coi_2018-20-v1/occurrence.txt", sep = "\t", header=TRUE) 
#ARMS_COI_2018_2020_eMOF <- read.csv("dwca-arms_coi_2018-20-v1/extendedmeasurementorfact.txt", sep = "\t", header=TRUE) 

#read taxonomic assignment files
April2021_COI_noBlank <- read.csv("tax_assignments_April2021_COI_noBlank.tsv", sep = "\t", header=FALSE) 
August2023_COI <- read.csv("tax_assignments_August2023_COI.tsv", sep = "\t", header=FALSE) 
COI <- read.csv("tax_assignments_COI.tsv", sep = "\t", header=FALSE) 
Jan2020_COI_noBlank <- read.csv("tax_assignments_Jan2020_COI_noBlank.tsv", sep = "\t", header=FALSE) 
Jan2022_COI_noBlank <- read.csv("tax_assignments_Jan2022_COI_noBlank.tsv", sep = "\t", header=FALSE) 
July2019_COI <- read.csv("tax_assignments_July2019_COI.tsv", sep = "\t", header=FALSE) 
May2021_COI <- read.csv("tax_assignments_May2021_COI.tsv", sep = "\t", header=FALSE) 
Sept2020_COI_noBlank <- read.csv("tax_assignments_Sept2020_COI_noBlank.tsv", sep = "\t", header=FALSE) 
#merge taxonomic assignment files
tax_assignments <- rbind(April2021_COI_noBlank, August2023_COI, COI, 
                         Jan2020_COI_noBlank, Jan2022_COI_noBlank, July2019_COI, 
                         May2021_COI, Sept2020_COI_noBlank)
#delete taxonomic assignment files
rm(April2021_COI_noBlank, August2023_COI, COI, 
  Jan2020_COI_noBlank, Jan2022_COI_noBlank, July2019_COI, 
  May2021_COI, Sept2020_COI_noBlank)
#delete columns with only NA values
tax_assignments <- tax_assignments[, colSums(is.na(tax_assignments)) != nrow(tax_assignments)]
#split V1 column into ASVID and abundance
tax_assignments[c('ASVID', 'abundance')] <- str_split_fixed(tax_assignments$V1, '_', 2)
#delete abundance column and the V12 column
tax_assignments <- select(tax_assignments, -abundance, -V12)
#correct the _ in the species names
tax_assignments$V20 <- gsub('_', ' ', tax_assignments$V20)
#rename the columns
colnames <- c("ID", "Kingdom", "Kingdom_bv", "Phylum", "Phylum_bv", "Class", "Class_bv", 
              "Order", "Order_bv", "Family", "Family_bv", "Genus", "Genus_bv", 
              "Species", "Species_bv", 'ASVID')
colnames(tax_assignments) <- colnames
rm(colnames)
# replace empty cells with NA
tax_assignments[tax_assignments == ''] <- NA

#delete the ID column
#tax_assignments_less <- select(tax_assignments, -ID)
#select the unique entries
#tax_assignments_less <- unique(tax_assignments_less)
#convert the ID column into the rowname of the data frame
#tax_assignments_less <- tax_assignments_less  %>% remove_rownames %>% column_to_rownames(var="ASVID")
#check the confidence estimates and if they are lower than 0.80 replace with NA
#this function will change the rownames, therefore we need to keep them in a vector and restore them later
rownames <- rownames(tax_assignments) 
tax_assignments <- replace_with_na_all(data = tax_assignments,
                                             condition = ~.x < 0.8)
tax_assignments <- as.data.frame(tax_assignments)
row.names(tax_assignments) <- rownames  

#correct the tax_assignments
taxonomy <- tax_assignments

#correct Phylum column
for (i in 1:nrow(taxonomy))
  if (is.na(taxonomy$Phylum_bv[i])){
    taxonomy$Phylum[i]=NA
  }

#correct Class column
for (i in 1:nrow(taxonomy))
  if (is.na(taxonomy$Class_bv[i])){
    taxonomy$Class[i]=NA
  }


#correct Order column
for (i in 1:nrow(taxonomy))
  if (is.na(taxonomy$Order_bv[i])){
    taxonomy$Order[i]=NA
  }


#correct Family column
for (i in 1:nrow(taxonomy))
  if (is.na(taxonomy$Family_bv[i])){
    taxonomy$Family[i]=NA
  }


#correct Genus column 
for (i in 1:nrow(taxonomy))
  if (is.na(taxonomy$Genus_bv[i])){
    taxonomy$Genus[i]=NA
  }


#correct Species column
for (i in 1:nrow(taxonomy))
  if (is.na(taxonomy$Species_bv[i])){
    taxonomy$Species[i]=NA
  }


tax_assignments <- taxonomy
rm(taxonomy)
rm(rownames, i)

#delete the columns with the confidence estimates
tax_assignments <- select(tax_assignments, -Kingdom_bv, -Phylum_bv, -Class_bv,
                          -Order_bv, -Family_bv, -Genus_bv, -Species_bv)

#correct the Kingdom column
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Phylum[i])){
    tax_assignments$Kingdom[i]='Biota incertae sedis'
  }

#add a column scientificName to keep the lowest taxonomic identification and the taxonRank
tax_assignments <- tax_assignments %>%
  mutate(scientificName = NA, 
         taxonRank = NA)
#populate those columns
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Species[i])){
    tax_assignments$scientificName[i]=tax_assignments$Genus[i]
    tax_assignments$taxonRank[i]='Genus'
  } else {
    tax_assignments$scientificName[i]=tax_assignments$Species[i]
    tax_assignments$taxonRank[i]='Species'
  }
    
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Genus[i])){
    tax_assignments$scientificName[i]=tax_assignments$Family[i]
    tax_assignments$taxonRank[i]='Family'
  }
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Family[i])){
    tax_assignments$scientificName[i]=tax_assignments$Order[i]
    tax_assignments$taxonRank[i]='Order'
  }
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Order[i])){
    tax_assignments$scientificName[i]=tax_assignments$Class[i]
    tax_assignments$taxonRank[i]='Class'
  }
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Class[i])){
    tax_assignments$scientificName[i]=tax_assignments$Phylum[i]
    tax_assignments$taxonRank[i]='Phylum'
  }
for (i in 1:nrow(tax_assignments))
  if (is.na(tax_assignments$Phylum[i])){
    tax_assignments$scientificName[i]=tax_assignments$Kingdom[i]
    tax_assignments$taxonRank[i]='Kingdom'
  }

#delete the detailed taxonomic levels
tax_assignments <- select(tax_assignments, -Kingdom, -Phylum, -Class, 
                          -Order, -Family, -Genus, -Species)

#correct the _ in the scientificName column
tax_assignments$scientificName <- gsub('_', ' ', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Hematodinium sp', 'Hematodinium', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Cystobasidium sp', 'Cystobasidium', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Cladosporium sp', 'Cladosporium', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Cladosporiumhaerospermum', 'Cladosporium sphaerospermum', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Sterigmatomyces sp', 'Sterigmatomyces', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Acanthamoeba sp', 'Acanthamoeba', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ectocarpus sp', 'Ectocarpus', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Lithothamnion sp', 'Lithothamnion', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Tolypocladium sp', 'Tolypocladium', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ilyonectria sp', 'Ilyonectria', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Eulalia aurea Gravier 1896', 'Eulalia aurea', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ischyrocerus anguipes 2', 'Ischyrocerus anguipes', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Metschnikowia 2', 'Metschnikowia', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ophelina acuminata CMC01', 'Ophelina acuminata', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ophelina acuminata CMC02', 'Ophelina acuminata', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Scalibregma inflatum CMC03', 'Scalibregma inflatum', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Apionsoma misakianum', 'Apionsoma (Apionsoma) misakianum', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Callyspongia siphonella', 'Callyspongia (Callyspongia) siphonella', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ciona roulei', 'Ciona roulii', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Clathria toxitenuis', 'Clathria (Microciona) toxitenuis', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Haliclona fascigera', 'Haliclona (Reniera) fascigera', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Halisarca dujardini', 'Halisarca dujardinii', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Myriotrichia claviformis', 'Myriotrichia clavaeformis', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Phyllophora pseudoceranoides', 'Phyllophora pseudoceranoïdes', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Polyides rotundus', 'Polyides rotunda', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Strongylacidon bermudae', 'Strongylacidon bermuda', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Veneroida', 'Venerida', tax_assignments$scientificName)
tax_assignments$scientificName <- gsub('Ocenebra erinacea', 'Ocinebra erinacea', tax_assignments$scientificName)

#select only the scientificName column for the WoRMS taxon match
scientificName <- select(tax_assignments, scientificName)
#select the unique ones
scientificName <- unique(scientificName)

#create a list
scientificName_list <- list()
for (i in scientificName) {
  scientificName_list <- append(scientificName_list, i)
}
scientificName_list <- as.character(scientificName_list)

#RETRIEVE APHIA IDs
#be careful here, run just the next line of code
scientificName_Aphia <- match_taxa(scientificName_list, ask=TRUE)
#the function will ask you if it finds ambihuous names
#you will need resolve them by typing the correct numbers 
#corresponding to the correct species names

#2  17968 Bdelloidea     Hudson, 1884 accepted exact     
#1  144562 Ceramium secundatum Lyngbye, 1819            accepted   exact     
#3  139129 Chiton olivaceus              L. Spengler, 1797 unaccepted             exact         
#1  167698 Clathria (Microciona) toxitenuis   Topsent, 1925               accepted   exact_subgenus
#2  394301 Clathria prolifera            Burton, 1940             unaccepted exact         
#1  137916 Doto           Oken, 1815    accepted       exact     
#1  145404 Ectocarpus fasciculatus Harvey, 1841  accepted   exact     
#1  129296 Glycera        Lamarck, 1818 accepted   exact     
#1  126892 Gobius niger   Linnaeus, 1758 accepted   exact     
#1  156508 Haliclona elegans (Lendenfeld, 1887) accepted   exact     
#2  109540 Heterocapsa    F.Stein, 1883 accepted exact     
#1  102431 Jassa falcata  (Montagu, 1808)                             accepted       exact     
#1  100182 Metschnikowia  T. Kamienski, 1899 accepted exact     
#1  121104 Monoposthia costata (Bastian, 1865) de Man, 1889                     accepted   exact     
#1  149142 Navicula       J.B.M. Bory de Saint-Vincent, 1822 accepted                 exact     
#1  149045 Nitzschia      A.H. Hassall, 1845 accepted exact     
#2  379195 Nymphon gracile Hoek, 1881          unaccepted exact     
#1  117034 Obelia         Péron & Lesueur, 1810 accepted   exact     
#1  140744 Philine aperta (Linnaeus, 1767)                             accepted   exact     
#1  143853 Polysiphonia   Greville, 1823       accepted   exact     
#2  291409 Pseudonitzschia NA             unaccepted exact     
#1  109572 Symbiodinium   LaJeunesse, 2017              accepted   exact     
#1     982 Terebellidae   Johnston, 1846            accepted   exact     
#1  100245 Tetracladium   De Wild., 1893   accepted exact     
#1 1870349 Trichoderma    Pers. 1794   accepted   exact     

#delete the rows for which no aphiaIDs were found
scientificName_Aphia <- scientificName_Aphia %>% drop_na(scientificName)

#delete match_type column and acceptedNameUsageID column
scientificName_Aphia <- select(scientificName_Aphia, -match_type, -acceptedNameUsageID)

#check which taxa don't have aphiaIDs
non_aphia <- anti_join(scientificName, scientificName_Aphia) 
#add scientificNameID column
non_aphia <- non_aphia %>%
  mutate(scientificNameID = 'no match')
#merge the two tables
scientificName_Aphia <- rbind(scientificName_Aphia, non_aphia)
rm(non_aphia)
rm(scientificName)

#merge it with the tax_assignments table
tax_assignments <- merge(scientificName_Aphia, tax_assignments)
#check if there are empty columns in tax_assignments
colSums(is.na(tax_assignments))
rm(scientificName_Aphia)
rm(i, scientificName_list)


################################################################################
################################################################################
################################################################################

#split taxonID column into ASVnumber and ASVID
ARMS_COI_2018_2020_occurrence[c('ASVnumber', 'ASVID')] <- str_split_fixed(ARMS_COI_2018_2020_occurrence$taxonID, ':', 2)
#delete ASVnumber column
ARMS_COI_2018_2020_occurrence <- select(ARMS_COI_2018_2020_occurrence, -ASVnumber)
#delete the columns that will be replaced
ARMS_COI_2018_2020_occurrence <- select(ARMS_COI_2018_2020_occurrence, -scientificName, -taxonRank, 
                                        -scientificNameID)
#delete unwanted columns
tax_assignments <- select(tax_assignments, -ID)

#select the ASVIDs from the tax_assignments that exist in the ARMS_COI_2018_2020_occurrence
matched_ASVIDs <- semi_join(tax_assignments, ARMS_COI_2018_2020_occurrence) 
#select the unique entries
matched_ASVIDs <- unique(matched_ASVIDs)

#merge the two files
ARMS_COI_2018_2020_occurrence <- merge(ARMS_COI_2018_2020_occurrence, matched_ASVIDs)
#delete the technical column
ARMS_COI_2018_2020_occurrence <- select(ARMS_COI_2018_2020_occurrence, -ASVID)

write.table(ARMS_COI_2018_2020_occurrence, "dwca-arms_coi_2018-20-v1/occurrence_corrected.txt", quote = FALSE,na = "",
            row.names = FALSE, col.names = TRUE, sep = "\t")

################################################################################
################################################################################
################################################################################

save.image("ARMS_issues_COI_2018_2020.RData")
#creating ".RData" in current working directory
#to save the basic files

#Now you have everything in your computer, and you can load it anytime you want by running
#load("ARMS_issues_COI_2018_2020.RData")
