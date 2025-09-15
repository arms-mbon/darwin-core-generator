"""
TODO series = df.iloc[0]
series are indexable like a dict
"""
import json
import traceback
import arms_toolbox.function as fn
from darwin_core.dwca import DwCArchive
from darwin_core.model import Table
from tqdm import tqdm
from .remote_data_file import ObservatoryDataFrame, OmicsDataFrame, SamplingDataFrame, MDAFasta, PEMADataFrame
import logging

logger = logging.getLogger(__name__)

class Pipeline:
    def __init__(self, genomic_region, time_window, aligned_assignment_url, data_release):
        # parameters
        self.genomic_region = genomic_region
        self.time_window = time_window
        self.aligned_assignment_url = aligned_assignment_url
        self.data_release = data_release

        # input
        self.observatory = ObservatoryDataFrame().fetch() # TODO observatory_superset
        self.omics = OmicsDataFrame().fetch()
        self.sampling = SamplingDataFrame().fetch()
        self.pema = PEMADataFrame(self.genomic_region, self.time_window, self.data_release).fetch()
        self.mda = MDAFasta(self.genomic_region, self.aligned_assignment_url).fetch()
        
        # output
        self.occurrence_core = Table(schema_path=f"./data/schemas/occurrence_core_schema_{self.genomic_region.lower()}.json")
        self.dna_extension = Table(schema_path=f"./data/schemas/dna_extension_schema_{self.genomic_region.lower()}.json")
        self.extended_measurement_or_facts_extension = Table(schema_path="./data/schemas/extended_measurement_or_facts_extension_schema.json")
        self.emof_properties = json.load(open("./data/schemas/extended_measurement_or_facts_extension_properties_schema.json"))
        self.emof_functions = {
            "submergedTime": fn.submerged_time,
            "preservativeUsed": fn.preservative_used,
            "lowerLimitFilterSize": fn.lower_limit_filter_size,
            "fieldReplicateID": fn.field_replicate,
            "technicalReplicateID": fn.technical_replicate,
        }
    
    def run(self):
        logger.info(f"\t->\trunning Pipeline{self.genomic_region}|{self.time_window}|{self.aligned_assignment_url}|{self.data_release}")
        self.replicate_material_sample_ids = sorted(list(set(self.pema.columns) - {"OTU", "Classification", "TAXON:NCBI_TAX_ID"}))
        self._set_context()
        if self.genomic_region in ("ITS", "18S"):
            self._map_specific_its_18s()
        elif self.genomic_region == "COI":
            self._map_specific_coi()
        else:
            raise AssertionError
        self._write()

    def _set_context(self):

        # slice DFs
        ...

        # web apis should also go here?
        ...
    
    

    def _map_specific_its_18s(self):
        for rmsid in tqdm(self.replicate_material_sample_ids):
            df_sampling_subset = self.sampling[self.sampling["ReplicateMaterialSampleID"] == rmsid]
            if len(df_sampling_subset) == 0:
                print(f"{rmsid} not found in sampling event data")
                continue
            if len(df_sampling_subset) > 1:
                print(f"{rmsid} is duplicated in sampling event data") # in this case, this code takes the first one only
                continue

            df_omics_subset = self.omics[self.omics["ReplicateMaterialSampleID"] == rmsid]
            if len(df_omics_subset) == 0:
                print(f"{rmsid} not found in omics data")
                continue
            if len(df_omics_subset) > 1:
                print(f"{rmsid} is duplicated in omics data") # in this case, this code takes the first one only
                continue

            oid = df_sampling_subset["ObservatoryID"].iloc[0]
            uid = df_sampling_subset["UnitID"].iloc[0]
            df_observatory_subset = self.observatory[(self.observatory["ObservatoryID"] == oid) & (self.observatory["UnitID"] == uid)]
            if len(df_observatory_subset) == 0:
                print(f"{oid}_{uid} not found in observatory data")
                continue
            if len(df_observatory_subset) > 1:
                print(f"{oid}_{uid} is duplicated in observatory data")
                continue

            df_pema_subset = self.pema[self.pema[rmsid] > 1]
            for _, row in df_pema_subset.iterrows():
                node_id = f"{rmsid}:{row['OTU']}"
                node_label = row["TAXON:NCBI_TAX_ID"]
                try:
                    self.occurrence_core.add_row(
                        occurrenceID=node_id,
                        eventID=df_sampling_subset["EventID"].iloc[0],
                        materialSampleID = fn.material_sample_id(df_omics_subset[f"Gene_{self.genomic_region}"].iloc[0]),
                        eventDate=f"{df_sampling_subset['DateDeployed'].iloc[0]}/{df_sampling_subset['DateCollected'].iloc[0]}",
                        organismQuantity=row[rmsid],
                        sampleSizeValue=df_pema_subset[rmsid].sum(),
                        associatedSequences=df_omics_subset[f"Gene_{self.genomic_region}"].iloc[0],
                        locationID=df_observatory_subset["MarineRegion_larger"].iloc[0],
                        decimalLatitude=df_observatory_subset["Latitude"].iloc[0],
                        decimalLongitude=df_observatory_subset["Longitude"].iloc[0],
                        countryCode=fn.country_code(df_observatory_subset),
                        fieldNumber=f"{oid}_{uid}",
                        minimumDepthInMeters=df_observatory_subset["Depth_min"].iloc[0],
                        maximumDepthInMeters=df_observatory_subset["Depth_max"].iloc[0],
                        habitat=fn.habitat(df_observatory_subset),
                        occurrenceRemarks=fn.occurrence_remarks(node_label),
                        eventRemarks=fn.event_remarks(df_sampling_subset),
                        verbatimIdentification=fn.verbatim_identification(row["Classification"]),
                        scientificNameID=fn.scientific_name_id(node_label, row["Classification"]),
                        scientificName=fn.scientific_name(node_label, row["Classification"]),
                        taxonRank=fn.taxon_rank(node_label, row["Classification"]),
                        taxonId=row["OTU"],
                        taxonConcept=f"NCBI:{node_label.split(':')[-1]}",
                        catalogNumber=df_sampling_subset["ReplicateMaterialSampleID"].iloc[0],
                    )

                    self.dna_extension.add_row(
                        occurrenceID=node_id,
                        env_broad_scale=df_observatory_subset["ENVO broad scale"].iloc[0],
                        env_local_scale=df_observatory_subset["ENVO local scale"].iloc[0],
                        env_medium_scale=df_observatory_subset["ENVO medium scale"].iloc[0],
                        DNA_Sequence=fn.dna_sequence(row["OTU"], self.mda),
                    )

                    for p in self.emof_properties:
                        f = self.emof_functions.get(p["measurementType"])
                        if f:
                            mv, mvid = f(
                                taxon_ncbi_tax_id=node_label,
                                df_sampling_subset=df_sampling_subset,
                                df_omics_subset=df_omics_subset,
                            )
                        else:
                            mv = p["measurementValue"]
                            mvid = p["measurementValueID"]
                        self.extended_measurement_or_facts_extension.add_row(
                            occurrenceID=node_id,
                            measurementType=p["measurementType"],
                            measurementUnit=p["measurementUnit"],
                            measurementValue=mv,
                            measurementTypeID=p["measurementTypeID"],
                            measurementUnitID=p["measurementUnitID"],
                            measurementValueID=mvid,
                        )
                except Exception as e:
                    print(f"{node_id} failed with exception {e}")
                    traceback.print_exc()

        
    def _map_specific_coi(self):
        for rmsid in tqdm(self.replicate_material_sample_ids):
            df_sampling_subset = self.sampling[self.sampling["ReplicateMaterialSampleID"] == rmsid]
            if len(df_sampling_subset) == 0:
                print(f"{rmsid} not found in sampling event data")
                continue
            if len(df_sampling_subset) > 1:
                print(f"{rmsid} is duplicated in sampling event data") # in this case, this code takes the first one only
                continue

            df_omics_subset = self.omics[self.omics["ReplicateMaterialSampleID"] == rmsid]
            if len(df_omics_subset) == 0:
                print(f"{rmsid} not found in omics data")
                continue
            if len(df_omics_subset) > 1:
                print(f"{rmsid} is duplicated in omics data") # in this case, this code takes the first one only
                continue

            oid = df_sampling_subset["ObservatoryID"].iloc[0]
            uid = df_sampling_subset["UnitID"].iloc[0]
            df_observatory_subset = self.observatory[(self.observatory["ObservatoryID"] == oid) & (self.observatory["UnitID"] == uid)]
            if len(df_observatory_subset) == 0:
                print(f"{oid}_{uid} not found in observatory data")
                continue
            if len(df_observatory_subset) > 1:
                print(f"{oid}_{uid} is duplicated in observatory data")
                continue

            df_pema_subset = self.pema[self.pema[rmsid] > 1]
            for _, row in df_pema_subset.iterrows():
                node_id = f"{rmsid}:{row['OTU']}"
                node_label = row["TAXON:NCBI_TAX_ID"]
                try:
                    self.occurrence_core.add_row(
                        occurrenceID=node_id,
                        eventID=df_sampling_subset["EventID"].iloc[0],
                        materialSampleID = fn.material_sample_id(df_omics_subset[f"Gene_{self.genomic_region}"].iloc[0]),
                        eventDate=f"{df_sampling_subset['DateDeployed'].iloc[0]}/{df_sampling_subset['DateCollected'].iloc[0]}",
                        organismQuantity=row[rmsid],
                        sampleSizeValue=df_pema_subset[rmsid].sum(),
                        associatedSequences=df_omics_subset[f"Gene_{self.genomic_region}"].iloc[0],
                        locationID=df_observatory_subset["MarineRegion_larger"].iloc[0],
                        decimalLatitude=df_observatory_subset["Latitude"].iloc[0],
                        decimalLongitude=df_observatory_subset["Longitude"].iloc[0],
                        countryCode=fn.country_code(df_observatory_subset),
                        fieldNumber=f"{oid}_{uid}",
                        minimumDepthInMeters=df_observatory_subset["Depth_min"].iloc[0],
                        maximumDepthInMeters=df_observatory_subset["Depth_max"].iloc[0],
                        habitat=fn.habitat(df_observatory_subset),
                        occurrenceRemarks=fn.occurrence_remarks(node_label),
                        eventRemarks=fn.event_remarks(df_sampling_subset),
                        verbatimIdentification=fn.verbatim_identification(row["Classification"]),
                        scientificNameID=fn.scientific_name_id(node_label, row["Classification"]),
                        scientificName=fn.scientific_name(node_label, row["Classification"]),
                        taxonRank=fn.taxon_rank(node_label, row["Classification"]),
                        taxonId=row["OTU"],
                        taxonConcept=f"NCBI:{node_label.split(':')[-1]}",
                        catalogNumber=df_sampling_subset["ReplicateMaterialSampleID"].iloc[0],
                    )

                    self.dna_extension.add_row(
                        occurrenceID=node_id,
                        env_broad_scale=df_observatory_subset["ENVO broad scale"].iloc[0],
                        env_local_scale=df_observatory_subset["ENVO local scale"].iloc[0],
                        env_medium_scale=df_observatory_subset["ENVO medium scale"].iloc[0],
                        DNA_Sequence=fn.dna_sequence(row["OTU"], self.mda, grouped=True),
                    )

                    for p in self.emof_properties:
                        f = self.emof_functions.get(p["measurementType"])
                        if f:
                            mv, mvid = f(
                                row=row,
                                taxon_ncbi_tax_id=node_label,
                                df_sampling_subset=df_sampling_subset,
                                df_omics_subset=df_omics_subset,
                                ptax=self.ptax # COI-specific, will probably return None for ITS/18S
                            )
                        else:
                            mv = p["measurementValue"]
                            mvid = p["measurementValueID"]
                        self.extended_measurement_or_facts_extension.add_row(
                            occurrenceID=node_id,
                            measurementType=p["measurementType"],
                            measurementUnit=p["measurementUnit"],
                            measurementValue=mv,
                            measurementTypeID=p["measurementTypeID"],
                            measurementUnitID=p["measurementUnitID"],
                            measurementValueID=mvid,
                        )
                except Exception as e:
                    print(f"{node_id} failed with exception {e}")
                    traceback.print_exc()

    def _write(self):
        dwca = DwCArchive()
        dwca.add_core(self.occurrence_core, "occurrence.csv")
        dwca.add_extension(self.dna_extension, "dnaextension.csv")
        dwca.add_extension(self.extended_measurement_or_facts_extension, "emof.csv")
        dwca.write(f"./data/output/{self.time_window}_{self.genomic_region}")
