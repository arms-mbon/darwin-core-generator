"""
Read in the ARMS-MBON data input -- PEMA files (OTU tables), setup files (for the emof), and observatory, omics, and sampling event 
 files -- and turn those into DwCA (occurrence, emof, and dnaextension CSV files).

To be improved
- Having the NCBI Id - APhia ID search put results in a dictionary so that they don't need to be done each time (? unsure if already done?)
- Reading from the tax_assigments files for COI: all our results go to species level, but this may not be so need to take the scientific name, 
  confidence level, and rank, from the entry that is the most specific for that row (not assume that it will be always the last row)
- There are now so many IF COI ELSE in here, that it needs to be split into two scripts (as was originally the case)
- I get reports for 18S that there are run accession numbers with no sample accession numbers found in ENA
  in fact that is not the case, it misses only some of the occurrence IDs
  """

import sys
import json
import pandas as pd
import requests
import time
import traceback
from Bio import SeqIO
from datetime import datetime
from io import StringIO
from pycountry import countries
from darwin_core.dwca import DwCArchive
from darwin_core.model import OccurrenceCore, DNAExtension, ExtendedMeasurementOrFactsExtension
from .pipeline import Pipeline

class PipelineCOI(Pipeline):
    def __init__(self, time_window, aligned_assignment_url):
        self.time_window = time_window
        self.aligned_assignment_url = aligned_assignment_url.replace("directlink.php?fid", "download.php?file")
        self.genomic_region = "COI"
        self.occurrence_core_schema = "./data/schemas/occurrence_core_schema_coi.json"
        self.dna_extension_schema = "./data/schemas/dna_extension_schema_coi.json"
        self.cache = {}
        self.ncbi2aphia = {}
        self.ncbi_multiples = {}  # occurrence_ids with multiple ncbi_tax_ids
        self.aphia_multiples = {}  # ncbi_tax_ids with multiple aphia_ids
    
    def run(self):
        # Read in the PEMA files from GH data_workspace repo
        full_file_name = f"https://raw.githubusercontent.com/arms-mbon/data_workspace/main/analysis_data/from_pema/processing_batch1/updated_taxonomic_assignments/Extended_final_table_{self.time_window}_{self.genomic_region}_noBlank_TaxonomyFull.csv"
        df_pema = pd.read_csv(full_file_name).rename(columns={"ASV_number:amplicon": "OTU"})            

        # Reading in the taxonomic assigments file, from where I want to get the confidence level for the final name 
        df_pema_tax_url = f"https://raw.githubusercontent.com/arms-mbon/data_workspace/main/analysis_data/from_pema/processing_batch1/taxonomic_assignments/tax_assignments_{self.time_window}_{self.genomic_region}_noBlank.tsv"
        df_pema_tax = pd.read_csv(df_pema_tax_url, sep=r"\s+", header=None)
        df_pema_tax[0] = df_pema_tax[0].apply(lambda _: _.split("_")[0]) # the first col are IDs, of which we want only the first part, to match to those IDs in the extended table
        df_pema_tax.columns = [
            "amplicon",
            "kingdom", "kingdom_confidence",
            "phylum", "phylum_confidence",
            "class", "class_confidence",
            "order", "order_confidence",
            "g",
            "family", "family_confidence",
            "genus", "genus_confidence",
            "species", "species_confidence",
        ]

        # Read in the observatory, omics, and sample event files from GH
        df_observatory = pd.read_csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_ObservatoryData.csv")
        df_omics = pd.read_csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_OmicsData.csv")
        df_sampling = pd.read_csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_SamplingEventData.csv")
        replicate_material_sample_ids = list(set(df_pema.columns) - {"OTU", "Classification", "TAXON:NCBI_TAX_ID"})
        replicate_material_sample_ids.sort()

        # Set up to read in fasta files (which are on the MDA so need to be downloaded)
        fasta = self.parse_fasta(self.aligned_assignment_url)

        # To get specific details on COI species results
        # For our data we have always species - later modify to look for most specific one available
        def get_species_confidence_level(otuid):
            try:
                matchrow = df_pema_tax.loc[df_pema_tax["amplicon"] == otuid]
                species_name = matchrow["species"].values[0]
                conflevel = matchrow["species_confidence"].values[0]
                return ( species_name, conflevel,"species")
            except Exception as E:
                print("missing match in the tax_assigments for ",otuid)
                return ("none","0","none")


        # Set up the rows for each occurrence that will be added to the emof
        emof_properties = json.load(open("./data/schemas/extended_measurement_or_facts_extension_properties_schema_coi.json"))
        emof_functions = {
            "submergedTime": self.get_emof_submerged_time,
            "preservativeUsed": self.get_emof_preservative_used,
            "lowerLimitFilterSize": self.get_emof_lower_limit_filter_size,
            "sampleAccessionNumber": self.get_emof_sample_accession_number,
            "NCBIID": self.get_emof_ncbi_id,
            "NCBIScientificName": self.get_emof_ncbi_scientific_name,
            "NCBITaxonRank": self.get_emof_ncbi_taxon_rank,
            "fieldReplicateID": self.get_emof_field_replicate,
            "technicalReplicateID": self.get_emof_technical_replicate,
        }

        # Set up the other CSV files of the DwCA
        occurrence_core = OccurrenceCore(schema_path=self.occurrence_core_schema)
        dna_extension = DNAExtension(schema_path=self.dna_extension_schema)
        extended_measurement_or_facts_extension = ExtendedMeasurementOrFactsExtension(schema_path="./data/schemas/extended_measurement_or_facts_extension_schema.json")
        
        # Run thru the input file going sample by sample (column by columm of the input file), where the column name is the replicate material sample
        # Will only work on those values in the column that are >1 
        counter = 0
        print("  Iterting thru the samples")
        for rmsid in replicate_material_sample_ids:
            if counter % 10 == 0:
                print("   count: ",counter)
            counter+=1 
            # check that this sample exists in the observatory, omics, and event CSVs in GH
            df_sampling_subset = df_sampling[df_sampling["ReplicateMaterialSampleID"] == rmsid]
            if len(df_sampling_subset) == 0:
                print(f"{rmsid} not found in sampling event data")
                continue
            if len(df_sampling_subset) > 1:
                print(f"{rmsid} is duplicated in sampling event data") # in this case, this code takes the first one only
                continue

            df_omics_subset = df_omics[df_omics["ReplicateMaterialSampleID"] == rmsid]
            if len(df_omics_subset) == 0:
                print(f"{rmsid} not found in omics data")
                continue
            if len(df_omics_subset) > 1:
                print(f"{rmsid} is duplicated in omics data") # in this case, this code takes the first one only
                continue

            oid = df_sampling_subset["ObservatoryID"].iloc[0]
            uid = df_sampling_subset["UnitID"].iloc[0]
            df_observatory_subset = df_observatory[(df_observatory["ObservatoryID"] == oid) & (df_observatory["UnitID"] == uid)]
            if len(df_observatory_subset) == 0:
                #print(f"{oid}_{uid} not found in observatory data")
                continue
            if len(df_observatory_subset) > 1:
                #print(f"{oid}_{uid} is duplicated in observatory data")
                continue

            # For each value in the replicate sample column that is >1, need to extract its related information from the different source files and 
            # add that to the DwC files (emof, occurrence, dnaextension)
            df_pema_subset = df_pema[df_pema[rmsid] > 1]
            for _, row in df_pema_subset.iterrows():
                occurrence_id = f"{rmsid}:{row['OTU']}"
                taxon_ncbi_tax_id = row["TAXON:NCBI_TAX_ID"]
                asvid_from_extended = row["OTU"].split(":")[1]
                try:
                    self.update_cache(occurrence_id, taxon_ncbi_tax_id)
                    occurrence_core.add_row(
                        occurrenceID=occurrence_id,
                        eventID=df_sampling_subset["EventID"].iloc[0],
                        materialSampleID=df_sampling_subset["ReplicateMaterialSampleID"].iloc[0],
                        eventDate=f"{df_sampling_subset['DateDeployed'].iloc[0]}/{df_sampling_subset['DateCollected'].iloc[0]}",
                        organismQuantity=row[rmsid],
                        sampleSizeValue=df_pema_subset[rmsid].sum(),
                        associatedSequences=df_omics_subset[f"Gene_{self.genomic_region}"].iloc[0],
                        locationID=df_observatory_subset["MarineRegion_larger"].iloc[0],
                        decimalLatitude=df_observatory_subset["Latitude"].iloc[0],
                        decimalLongitude=df_observatory_subset["Longitude"].iloc[0],
                        countryCode=countries.get(alpha_3=df_observatory_subset["Country"].iloc[0]).alpha_2,
                        fieldNumber=f"{oid}_{uid}",
                        minimumDepthInMeters=df_observatory_subset["Depth_min"].iloc[0],
                        maximumDepthInMeters=df_observatory_subset["Depth_max"].iloc[0],
                        habitat=self.get_core_habitat(df_observatory_subset),
                        occurrenceRemarks=self.get_core_occurrence_remarks(taxon_ncbi_tax_id),
                        eventRemarks=self.get_core_event_remarks(df_sampling_subset),
                        verbatimIdentification=self.get_core_verbatim_identification(row["Classification"]),
                        scientificNameID=self.cache[taxon_ncbi_tax_id]["scientificNameID"],
                        scientificName=self.cache[taxon_ncbi_tax_id]["scientificName"],
                        taxonRank=self.cache[taxon_ncbi_tax_id]["taxonRank"],
                        taxonId=row["OTU"],
                    )

                    dna_extension.add_row(
                        occurrenceID=occurrence_id,
                        env_broad_scale=df_observatory_subset["ENVO broad scale"].iloc[0],
                        env_local_scale=df_observatory_subset["ENVO local scale"].iloc[0],
                        env_medium_scale=df_observatory_subset["ENVO medium scale"].iloc[0],
                        DNA_Sequence=self.get_dna_dna_sequence(row["OTU"], fasta),
                    )

                    # Bram may want to fix this hacky code
                    # Adding in the extra 3 entries for the COI emof
                    otuid = row["OTU"].split(":")[1] # the ID in the fasta file is the second part of the ASV ID from the extended final table
                    name, confidence_level, rank = get_species_confidence_level(otuid) #XXX
                    emof_functions["originalScientificName"] = fnwrapper(name)
                    emof_functions["originalTaxonRank"]= fnwrapper(rank)
                    emof_functions["originalScientificNameConfidenceLevel"] = fnwrapper(confidence_level)

                    for p in emof_properties:
                        fn = emof_functions.get(p["measurementType"])
                        if fn:
                            mv, mvid = fn(
                                taxon_ncbi_tax_id=taxon_ncbi_tax_id,
                                df_sampling_subset=df_sampling_subset,
                                df_omics_subset=df_omics_subset,
                            )
                        else:
                            mv = p["measurementValue"]
                            mvid = p["measurementValueID"]
                        extended_measurement_or_facts_extension.add_row(
                            occurrenceID=occurrence_id,
                            measurementType=p["measurementType"],
                            measurementUnit=p["measurementUnit"],
                            measurementValue=mv,
                            measurementTypeID=p["measurementTypeID"],
                            measurementUnitID=p["measurementUnitID"],
                            measurementValueID=mvid,
                        )
                except Exception as e:
                    print(f"{occurrence_id} failed with exception {e}")
                    traceback.print_exc()

        dwca = DwCArchive()
        dwca.add_core(occurrence_core, "occurrence.csv")
        dwca.add_extension(dna_extension, "dnaextension.csv")
        dwca.add_extension(extended_measurement_or_facts_extension, "emof.csv")
        print("writing out to ./data/output/",{self.time_window},"_",{self.genomic_region}) # XXX
        print("")
        dwca.write(f"./data/output/{self.time_window}_{self.genomic_region}")
        self.report_multiples(
            ncbi_path=f"./data/output/{self.time_window}_{self.genomic_region}/ncbi_multiples.json",
            aphia_path=f"./data/output/{self.time_window}_{self.genomic_region}/aphia_multiples.json"
        )
        #if self.genomic_region == "ITS":  # hack to get rid of technicalReplicateID in ITS DECIDED TO KEEP
        #    df = pd.read_csv(f"./data/output/{self.time_window}_ITS/emof.csv")
        #    df = df[df["measurementType"] != "technicalReplicateID"]
        #    df.to_csv(f"./data/output/{self.time_window}_ITS/emof.csv", index=False)

    # Download the fasta file from MDA
    @staticmethod
    def parse_fasta(download_url):
        response = requests.get(download_url).content.decode()
        time.sleep(1)
        fasta_grouped = SeqIO.to_dict(SeqIO.parse(StringIO(response), "fasta"))
        fasta = {k.split("_")[0]: v for k, v in fasta_grouped.items()}
        assert len(fasta) == len(fasta_grouped)
        return fasta
    
    # get the taxon id (from the file column's title) being the second part of the string
    @staticmethod
    def parse_pema_taxon(taxon_ncbi_tax_id):
        taxon = taxon_ncbi_tax_id.split(":")[0]
        ncbi_tax_ids = taxon_ncbi_tax_id.split(":")[1].split("\\n")
        return taxon, ncbi_tax_ids
    
    # get the rank information for the species name from ENA
    @staticmethod
    def request_ena(ncbi_tax_id):
        response = requests.get(
            f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{ncbi_tax_id}?binomialOnly=false",
            headers={"Accept": "application/json"}
        )
        time.sleep(1)
        if response.status_code == 200:
            return response.json()
        else:
            return None

    # Get the Aphia ID for the NCBI ID
    def request_aphia(self, ncbi_tax_id):
        if not ncbi_tax_id in self.ncbi2aphia.keys():
            response = requests.get(
                f"https://www.marinespecies.org/rest/AphiaRecordByExternalID/{ncbi_tax_id}?type=ncbi",
                headers={"Accept": "application/json"}
            )
            time.sleep(1)
            if response.status_code == 200:
                self.ncbi2aphia[ncbi_tax_id] = response.json()
            else:
                self.ncbi2aphia[ncbi_tax_id] = {}

    def update_cache(self, occurrence_id, taxon_ncbi_tax_id):
        if not taxon_ncbi_tax_id in self.cache.keys():
            self.cache[taxon_ncbi_tax_id] = {}
            taxon, ncbi_tax_ids = self.parse_pema_taxon(taxon_ncbi_tax_id)
            self.cache[taxon_ncbi_tax_id]["NCBIScientificName"] = taxon
            self.cache[taxon_ncbi_tax_id]["ncbi_tax_ids"] = ncbi_tax_ids
            for ncbi_tax_id in ncbi_tax_ids:
                response_json = self.request_ena(ncbi_tax_id)
                if response_json:
                    self.cache[taxon_ncbi_tax_id]["NCBIID"] = ncbi_tax_id
                    self.cache[taxon_ncbi_tax_id]["NCBITaxonRank"] = response_json["rank"]
                    break
            assert self.cache[taxon_ncbi_tax_id]["NCBIID"]
            ncbi_tax_id = self.cache[taxon_ncbi_tax_id]["NCBIID"]

            self.request_aphia(ncbi_tax_id)
            self.cache[taxon_ncbi_tax_id]["AphiaID"] = self.ncbi2aphia[ncbi_tax_id].get("AphiaID", "no match")
            self.cache[taxon_ncbi_tax_id]["scientificNameID"] = self.ncbi2aphia[ncbi_tax_id].get("lsid", "no match")
            self.cache[taxon_ncbi_tax_id]["scientificName"] = self.ncbi2aphia[ncbi_tax_id].get("scientificname", "no match")
            self.cache[taxon_ncbi_tax_id]["taxonRank"] = self.ncbi2aphia[ncbi_tax_id].get("rank", "no match")
            self.cache[taxon_ncbi_tax_id]["isTerrestrial"] = self.ncbi2aphia[ncbi_tax_id].get("isTerrestrial", "no match")
        
        if len(self.cache[taxon_ncbi_tax_id]["ncbi_tax_ids"]) > 1:
            self.ncbi_multiples[occurrence_id] = {}
            for ncbi_tax_id in self.cache[taxon_ncbi_tax_id]["ncbi_tax_ids"]:
                self.request_aphia(ncbi_tax_id)
                self.ncbi_multiples[occurrence_id][ncbi_tax_id] = {}
                self.ncbi_multiples[occurrence_id][ncbi_tax_id]["AphiaID"] = self.ncbi2aphia[ncbi_tax_id].get("AphiaID", "N/A")
                self.ncbi_multiples[occurrence_id][ncbi_tax_id]["lsid"] = self.ncbi2aphia[ncbi_tax_id].get("lsid", "N/A")
                self.ncbi_multiples[occurrence_id][ncbi_tax_id]["scientificname"] = self.ncbi2aphia[ncbi_tax_id].get("scientificname", "N/A")
                self.ncbi_multiples[occurrence_id][ncbi_tax_id]["rank"] = self.ncbi2aphia[ncbi_tax_id].get("rank", "N/A")
                self.ncbi_multiples[occurrence_id][ncbi_tax_id]["isTerrestrial"] = self.ncbi2aphia[ncbi_tax_id].get("isTerrestrial", "N/A")

    @staticmethod
    def get_dna_dna_sequence(otu, fasta):
        return fasta[otu.split(":")[1]].seq

    def report_multiples(self, ncbi_path, aphia_path):
        with open(ncbi_path, "w") as f:
            json.dump(self.ncbi_multiples, f, indent=4)
        with open(aphia_path, "w") as f:
            json.dump(self.aphia_multiples, f, indent=4)

    @staticmethod
    def get_core_habitat(df_observatory_subset):
        return "; ".join([
            df_observatory_subset["Monitoring area"].iloc[0].strip(),
            df_observatory_subset["Anthropogenic influence"].iloc[0].strip(),
            df_observatory_subset["IUCN habitat type"].iloc[0].strip(),
        ])

    def get_core_occurrence_remarks(self, taxon_ncbi_tax_id):
        if str(self.cache[taxon_ncbi_tax_id]["isTerrestrial"]) == "1":
            return "uncertain classification (not marine)"
        else:
            return ""

    @staticmethod
    def get_core_event_remarks(df_sampling_subset):
        crate_cover = df_sampling_subset["CrateCover"].iloc[0]
        fraction = df_sampling_subset["Fraction"].iloc[0]
        return f"Was a crate cover used when retrieving the unit: {crate_cover}; Which component of the ARMS unit material was extracted: {fraction}"

    @staticmethod
    def get_core_verbatim_identification(classification):
        prefix = "Cellular organisms;"
        if classification.startswith(prefix):
            return classification[len(prefix):]
        else:
            return classification
    
    @staticmethod
    def get_emof_submerged_time(**kwargs):
        date_deployed = kwargs["df_sampling_subset"]["DateDeployed"].iloc[0]
        date_collected = kwargs["df_sampling_subset"]["DateCollected"].iloc[0]
        mv = (
            datetime.strptime(date_collected, "%Y-%m-%d")
            - datetime.strptime(date_deployed, "%Y-%m-%d")
        ).days
        return mv, "" 

    @staticmethod
    def get_emof_preservative_used(**kwargs):
        mv = kwargs["df_sampling_subset"]["Preservative"].iloc[0]
        mvid = {
            "etoh": "https://pubchem.ncbi.nlm.nih.gov/compound/702",
            "dmso": "https://pubchem.ncbi.nlm.nih.gov/compound/679"
        }.get(mv.lower(), "")
        return mv, mvid

    @staticmethod
    def get_emof_lower_limit_filter_size(**kwargs):
        mv = kwargs["df_sampling_subset"]["Filter"].iloc[0]
        return mv, "" 

    # getting the sample accession number for a specific run accession number
    # adding in the dictionary that was still a TODO. associated_sequences (run accession) is the key and sample (SAMEA) is the value
    def get_emof_sample_accession_number(self, **kwargs):
        dictofsequences = {} 
        try:
            associated_sequences = kwargs["df_omics_subset"][f"Gene_{self.genomic_region}"].iloc[0]  # TODO store intermediate results in a dict {<associated_sequences>: <sample_accession_number>}
            if associated_sequences in dictofsequences.keys():
                mv = dictofsequences[associated_sequences]
            else:
                response = requests.post(
                    "https://www.ebi.ac.uk/ena/portal/api/search",
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    data={
                        "result": "read_run",
                        "fields": "sample_accession",
                        "includeAccessionType": "run",
                        "includeAccessions": associated_sequences,
                        "format": "tsv",
                    }
                )
                mv = pd.read_csv(StringIO(response.text), sep="\t")["sample_accession"].iloc[0]
                dictofsequences[associated_sequences] = mv 
        except:
            print("No sample accession found for ",associated_sequences)
            mv = "None"
        time.sleep(1)
        return mv,"" 

    def get_emof_ncbi_id(self, **kwargs):
        mv = self.cache[kwargs["taxon_ncbi_tax_id"]]["NCBIID"]
        return mv, "" 

    def get_emof_ncbi_scientific_name(self, **kwargs):
        mv = self.cache[kwargs["taxon_ncbi_tax_id"]]["NCBIScientificName"]
        return mv, "" 

    def get_emof_ncbi_taxon_rank(self, **kwargs):
        mv = self.cache[kwargs["taxon_ncbi_tax_id"]]["NCBITaxonRank"]
        return mv, "" 

    @staticmethod
    def get_emof_field_replicate(**kwargs):
        mv = kwargs["df_sampling_subset"]["FieldReplicate"].iloc[0]
        return mv, "" 

    def get_emof_technical_replicate(self, **kwargs):
        # Will do this for ITS after all, so changing next line
        #if (not (self.genomic_region == "ITS")) and (kwargs["df_sampling_subset"]["SequencingRunRepeat"].iloc[0] == "second sequencing run (repeat)"):
        if (kwargs["df_sampling_subset"]["SequencingRunRepeat"].iloc[0] == "second sequencing run (repeat)"):
            mv = kwargs["df_sampling_subset"]["MaterialSampleID"].iloc[0]
            return mv, "" 
        else:
            return "", "" 
        
def fnwrapper(mv):
    mvid = "" # because for these 3 terms we know the meas value id is blank. HACK
    # empty wrapper function
    def fn(*args,**kwargs):
        return mv,mvid
    return fn  