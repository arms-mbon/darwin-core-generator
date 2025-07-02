"""
this pipeline might be discarded in the future when the pema analysis is rerun
with the latest pema version
"""
import contextlib
import json
import pandas as pd
import requests
import time
import traceback
from Bio import SeqIO
from copy import deepcopy
from io import StringIO
from pycountry import countries
from darwin_core.dwca import DwCArchive
from darwin_core.model import OccurrenceCore, DNAExtension, ExtendedMeasurementOrFactsExtension
from arms_toolbox.pipeline_its import PipelineITS

class PipelineCOI(PipelineITS):
    def __init__(self, time_window, aligned_assignment_url):
        super().__init__(time_window, aligned_assignment_url)
        self.genomic_region = "COI"
        self.occurrence_core_schema = "./data/schemas/occurrence_core_schema_coi.json"
        self.dna_extension_schema = "./data/schemas/dna_extension_schema_coi.json"
        self.name2id = {}
    
    def run(self):
        df_pema = (
            pd.read_excel(f"https://raw.githubusercontent.com/arms-mbon/data_workspace/main/analysis_data/from_pema/processing_batch1/taxonomic_assignments/Extended_final_table_{self.time_window}_{self.genomic_region}_noBlank.xlsx")
            .rename(columns={"ASV_number:amplicon": "OTU"})
        )
        df_pema_tax = pd.read_csv(f"https://raw.githubusercontent.com/arms-mbon/data_workspace/main/analysis_data/from_pema/processing_batch1/taxonomic_assignments/tax_assignments_{self.time_window}_{self.genomic_region}_noBlank.tsv", sep=r"\s+", header=None)
        df_pema_tax[0] = df_pema_tax[0].apply(lambda _: _.split("_")[0])
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
        df_observatory = pd.read_csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_ObservatoryData.csv")
        df_omics = pd.read_csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_OmicsData.csv")
        df_sampling = pd.read_csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_SamplingEventData.csv")
        replicate_material_sample_ids = list(set(df_pema.columns) - {"OTU", "Classification", "TAXON:NCBI_TAX_ID"})
        replicate_material_sample_ids.sort()

        fasta = self.parse_fasta(self.aligned_assignment_url)

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
            "originalScientificName": self.get_emof_original_scientific_name,
            "originalTaxonRank": self.get_emof_original_taxon_rank,
            "originalScientificNameConfidenceLevel": self.get_emof_original_scientific_name_confidence_level,
        }

        occurrence_core = OccurrenceCore(schema_path=self.occurrence_core_schema)
        dna_extension = DNAExtension(schema_path=self.dna_extension_schema)
        extended_measurement_or_facts_extension = ExtendedMeasurementOrFactsExtension(schema_path="./data/schemas/extended_measurement_or_facts_extension_schema.json")

        for rmsid in replicate_material_sample_ids:
            df_sampling_subset = df_sampling[df_sampling["ReplicateMaterialSampleID"] == rmsid]
            if len(df_sampling_subset) == 0:
                print(f"{rmsid} not found in sampling event data")
                continue
            if len(df_sampling_subset) > 1:
                print(f"{rmsid} is duplicated in sampling event data")
                continue

            df_omics_subset = df_omics[df_omics["ReplicateMaterialSampleID"] == rmsid]
            if len(df_omics_subset) == 0:
                print(f"{rmsid} not found in omics data")
                continue
            if len(df_omics_subset) > 1:
                print(f"{rmsid} is duplicated in omics data")
                continue

            oid = df_sampling_subset["ObservatoryID"].iloc[0]
            uid = df_sampling_subset["UnitID"].iloc[0]
            df_observatory_subset = df_observatory[(df_observatory["ObservatoryID"] == oid) & (df_observatory["UnitID"] == uid)]
            if len(df_observatory_subset) == 0:
                print(f"{oid}_{uid} not found in observatory data")
                continue
            if len(df_observatory_subset) > 1:
                print(f"{oid}_{uid} is duplicated in observatory data")
                continue

            df_pema_subset = df_pema[df_pema[rmsid] > 1]

            for _, row in df_pema_subset.iterrows():
                occurrence_id = f"{rmsid}:{row['OTU']}"
                # print(occurrence_id)  # XXX
                try:
                    self.update_cache(row["OTU"], df_pema_tax)
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
                        occurrenceRemarks=self.get_core_occurrence_remarks(row["OTU"]),
                        eventRemarks=self.get_core_event_remarks(df_sampling_subset),
                        verbatimIdentification=self.get_core_verbatim_identification(row["Classification"]),
                        scientificNameID=self.cache[row["OTU"]]["scientificNameID"],
                        scientificName=self.cache[row["OTU"]]["scientificName"],
                        taxonRank=self.cache[row["OTU"]]["taxonRank"],
                        taxonId=row["OTU"],
                    )
                    dna_extension.add_row(
                        occurrenceID=occurrence_id,
                        env_broad_scale=df_observatory_subset["ENVO broad scale"].iloc[0],
                        env_local_scale=df_observatory_subset["ENVO local scale"].iloc[0],
                        env_medium_scale=df_observatory_subset["ENVO medium scale"].iloc[0],
                        DNA_Sequence=self.get_dna_dna_sequence(row["OTU"], fasta),
                    )
                    for p in emof_properties:
                        fn = emof_functions.get(p["measurementType"])
                        if fn:
                            mv, mvid = fn(
                                row=row,
                                df_sampling_subset=df_sampling_subset,
                                df_omics_subset=df_omics_subset,
                                df_pema_tax=df_pema_tax,
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
        dwca.write(f"./data/output/{self.time_window}_{self.genomic_region}")
        self.report_multiples(
            ncbi_path=f"./data/output/{self.time_window}_{self.genomic_region}/ncbi_multiples.json",
            aphia_path=f"./data/output/{self.time_window}_{self.genomic_region}/aphia_multiples.json"
        )

    @staticmethod
    def parse_fasta(download_url):
        response = requests.get(download_url).content.decode()
        time.sleep(1)
        fasta_grouped = SeqIO.to_dict(SeqIO.parse(StringIO(response), "fasta"))
        fasta = {k.split("_")[0]: v for k, v in fasta_grouped.items()}
        assert len(fasta) == len(fasta_grouped)
        return fasta
    
    def get_df_pema_tax_subset(self, df_pema_tax, otu):
        amplicon = otu.split(":")[1]
        df_pema_tax_subset = df_pema_tax[df_pema_tax["amplicon"] == amplicon]
        if len(df_pema_tax_subset) == 0:
            raise AssertionError(f"{amplicon} not found in tax assignments")
        if len(df_pema_tax_subset) > 1:
            raise AssertionError(f"{amplicon} is duplicated in tax assignments")
        return df_pema_tax_subset

    def fetch_id_from_name(self, name):
        if not name in self.name2id.keys():
            with contextlib.redirect_stdout(StringIO()):
                from arms_toolbox.ncbi_taxonomist_wrapper import resolve 
                import sys
                sys_exit_copy = deepcopy(sys.exit)
                sys.exit = lambda *args, **kwargs: None
                ncbi_id = resolve(name)
                sys.exit = sys_exit_copy
                time.sleep(2)  # XXX
                self.name2id[name] = ncbi_id
        return self.name2id[name]

    def update_cache(self, otu, df_pema_tax):
        if not otu in self.cache.keys():
            ncbi_id = None
            ncbi_scientific_name = None
            ncbi_taxon_rank = None
            original_scientific_name = None
            original_taxon_rank = None
            original_scientific_name_confidence_level = None
            df_pema_tax_subset = self.get_df_pema_tax_subset(df_pema_tax, otu)
            ranks = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
            for _rank in ranks:
                _original_name = df_pema_tax_subset[_rank].iloc[0]
                _name = _original_name.replace("_", " ")
                if _name.endswith(" sp"):
                    _name = f"unclassified {_name[:-3]}"
                _id = self.fetch_id_from_name(_name)
                if _id is None:
                    continue
                ncbi_scientific_name = _name
                ncbi_id = _id
                ncbi_taxon_rank = self.request_ena(_id)["rank"]
                original_scientific_name = _original_name
                original_taxon_rank = _rank
                original_scientific_name_confidence_level = df_pema_tax_subset[f"{_rank}_confidence"].iloc[0]
                break

            self.cache[otu] = {}
            self.cache[otu]["NCBIID"] = ncbi_id
            self.cache[otu]["NCBIScientificName"] = ncbi_scientific_name 
            self.cache[otu]["NCBITaxonRank"] = ncbi_taxon_rank
            self.cache[otu]["OriginalScientificName"] = original_scientific_name 
            self.cache[otu]["OriginalTaxonRank"] = original_taxon_rank
            self.cache[otu]["OriginalScientificNameConfidenceLevel"] = original_scientific_name_confidence_level

            self.request_aphia(ncbi_id)
            self.cache[otu]["AphiaID"] = self.ncbi2aphia[ncbi_id].get("AphiaID", "no match")
            self.cache[otu]["scientificNameID"] = self.ncbi2aphia[ncbi_id].get("lsid", "no match")
            self.cache[otu]["scientificName"] = self.ncbi2aphia[ncbi_id].get("scientificname", "no match")
            self.cache[otu]["taxonRank"] = self.ncbi2aphia[ncbi_id].get("rank", "no match")
            self.cache[otu]["isTerrestrial"] = self.ncbi2aphia[ncbi_id].get("isTerrestrial", "no match")

    @staticmethod
    def get_dna_dna_sequence(otu, fasta):
        return fasta[otu.split(":")[1]].seq

    def get_emof_ncbi_id(self, **kwargs):
        mv = self.cache[kwargs["row"]["OTU"]]["NCBIID"]
        return mv, "" 

    def get_emof_ncbi_scientific_name(self, **kwargs):
        mv = self.cache[kwargs["row"]["OTU"]]["NCBIScientificName"]
        return mv, "" 

    def get_emof_ncbi_taxon_rank(self, **kwargs):
        mv = self.cache[kwargs["row"]["OTU"]]["NCBITaxonRank"]
        return mv, "" 

    def get_emof_original_scientific_name(self, **kwargs):
        mv = self.cache[kwargs["row"]["OTU"]]["OriginalScientificName"]
        return mv, "" 

    def get_emof_original_taxon_rank(self, **kwargs):
        mv = self.cache[kwargs["row"]["OTU"]]["OriginalTaxonRank"]
        return mv, "" 

    def get_emof_original_scientific_name_confidence_level(self, **kwargs):
        mv = self.cache[kwargs["row"]["OTU"]]["OriginalScientificNameConfidenceLevel"]
        return mv, "" 
