"""
TODO separate transform from load
"""
import traceback
import arms_toolbox.function as fn
from tqdm import tqdm
from .data_source import ObservatoryDataFrame, OmicsDataFrame, SamplingDataFrame, MDAFasta, PEMADataFrame
from .data_sink import OccurrenceCore, DNAExtension, ExtendedMeasurementOrFactsExtension
from .parameter import EMOFProperties
import logging

logger = logging.getLogger(__name__)

class Pipeline:
    def __init__(self, genomic_region, time_window, aligned_assignment_url, data_release):
        # parameters
        self.genomic_region = genomic_region
        self.time_window = time_window
        self.aligned_assignment_url = aligned_assignment_url
        self.data_release = data_release
        self.emof_properties = EMOFProperties  # TODO refactor this
        self.emof_functions = {  # TODO refactor this
            "submergedTime": fn.submerged_time,
            "preservativeUsed": fn.preservative_used,
            "lowerLimitFilterSize": fn.lower_limit_filter_size,
            "fieldReplicateID": fn.field_replicate,
            "technicalReplicateID": fn.technical_replicate,
        }

        # input
        self.observatory_df = ObservatoryDataFrame().fetch()
        self.omics_df = OmicsDataFrame().fetch()
        self.sampling_df = SamplingDataFrame().fetch()
        self.pema_df = PEMADataFrame(self.genomic_region, self.time_window, self.data_release).fetch()
        self.mda_fasta = MDAFasta(self.genomic_region, self.aligned_assignment_url).fetch()
        
        # output
        self.occurrence_core = OccurrenceCore(self.genomic_region, self.time_window)
        self.dna_extension = DNAExtension(self.genomic_region, self.time_window)
        self.extended_measurement_or_facts_extension = ExtendedMeasurementOrFactsExtension(self.genomic_region, self.time_window)

    def run(self):
        logger.info(f"\t->\trunning Pipeline{self.genomic_region}|{self.time_window}|{self.aligned_assignment_url}|{self.data_release}")
        self._etl()
        self._serialize()
    
    def _etl(self):
        self.replicate_material_sample_ids = sorted(list(set(self.pema_df.columns) - {"OTU", "Classification", "TAXON:NCBI_TAX_ID"}))
        for rmsi in tqdm(self.replicate_material_sample_ids):
            self.rmsi = rmsi
            try:
                self._extract()
            except Exception as e:
                logger.error(f"ExtractionError; {self.rmsi} failed with exception {e}")
                # traceback.print_exc() # TODO log instead of print
                continue
            for _, pema in self.pema_df_gt.iterrows():
                self.pema = pema
                self.occurrence_id = f"{self.rmsi}:{self.pema['OTU']}"
                self.pema_tree = self.pema["Classification"]
                self.pema_leaf = self.pema["TAXON:NCBI_TAX_ID"]
                try:
                    if self.pema_tree.lower() == "not found" and self.pema_leaf.lower() == "not found":
                        raise AssertionError(f"{self.occurrence_id} has missing taxonomic classification")
                    self._transform_and_load_occurrence_core()
                    self._transform_and_load_dna_extension()
                    self._transform_and_load_extended_measurement_or_facts_extension()
                except Exception as e:
                    logger.error(f"TransformError; {self.occurrence_id} failed with exception {e}")
                    # traceback.print_exc() # TODO log instead of print
    
    def _extract(self):
        # sampling
        self.sampling = self.sampling_df[self.sampling_df["ReplicateMaterialSampleID"] == self.rmsi]
        if len(self.sampling) == 0:
            raise AssertionError(f"{self.rmsi} not found in sampling event data")
        if len(self.sampling) > 1:
            raise AssertionError(f"{self.rmsi} is duplicated in sampling event data")
        self.sampling = self.sampling.iloc[0]

        # omics
        self.omics = self.omics_df[self.omics_df["ReplicateMaterialSampleID"] == self.rmsi]
        if len(self.omics) == 0:
            raise AssertionError(f"{self.rmsi} not found in omics data")
        if len(self.omics) > 1:
            raise AssertionError(f"{self.rmsi} is duplicated in omics data")
        self.omics = self.omics.iloc[0]

        # observatory
        oid = self.sampling["ObservatoryID"]
        uid = self.sampling["UnitID"]
        self.observatory = self.observatory_df[(self.observatory_df["ObservatoryID"] == oid) & (self.observatory_df["UnitID"] == uid)]
        if len(self.observatory) == 0:
            raise AssertionError(f"combination of {oid} & {uid} not found in observatory data")
        if len(self.observatory) > 1:
            raise AssertionError(f"combination of {oid} & {uid} is duplicated in observatory data")
        self.observatory = self.observatory.iloc[0]

        # pema_df_gt (the others are a pd.Series, but this will still be a pd.DataFrame)
        self.pema_df_gt = self.pema_df[self.pema_df[self.rmsi] > 1]


    def _transform_and_load_occurrence_core(self):
        self.occurrence_core.add_node(
            occurrenceID=self.occurrence_id,  # node identifier
            eventID=self.sampling["EventID"],
            materialSampleID = fn.material_sample_id(self.omics[f"Gene_{self.genomic_region}"]),
            eventDate=f"{self.sampling['DateDeployed']}/{self.sampling['DateCollected']}", # TODO fn
            organismQuantity=self.pema[self.rmsi],
            sampleSizeValue=self.pema_df_gt[self.rmsi].sum(), # TODO fn
            associatedSequences=self.omics[f"Gene_{self.genomic_region}"],
            locationID=self.observatory["MarineRegion_larger"],
            decimalLatitude=self.observatory["Latitude"],
            decimalLongitude=self.observatory["Longitude"],
            countryCode=fn.country_code(self.observatory),
            fieldNumber=f"{self.sampling["ObservatoryID"]}_{self.sampling["UnitID"]}", # TODO fn
            minimumDepthInMeters=self.observatory["Depth_min"],
            maximumDepthInMeters=self.observatory["Depth_max"],
            habitat=fn.habitat(self.observatory),
            occurrenceRemarks=fn.occurrence_remarks(self.pema_leaf),
            eventRemarks=fn.event_remarks(self.sampling),
            verbatimIdentification=fn.verbatim_identification(self.pema_tree),
            scientificNameID=fn.scientific_name_id(self.pema_leaf, self.pema_tree),
            scientificName=fn.scientific_name(self.pema_leaf, self.pema_tree),
            taxonRank=fn.taxon_rank(self.pema_leaf, self.pema_tree),
            taxonId=self.pema["OTU"],
            taxonConcept=f"NCBI:{self.pema_leaf.split(':')[-1]}", # TODO fn
            catalogNumber=self.sampling["ReplicateMaterialSampleID"],
        )
    
    def _transform_and_load_dna_extension(self):
        self.dna_extension.add_node(
            occurrenceID=self.occurrence_id,  # node identifier
            env_broad_scale=self.observatory["ENVO broad scale"],
            env_local_scale=self.observatory["ENVO local scale"],
            env_medium_scale=self.observatory["ENVO medium scale"],
            DNA_Sequence=fn.dna_sequence(self.pema["OTU"], self.mda_fasta, self.genomic_region),
        )
    
    def _transform_and_load_extended_measurement_or_facts_extension(self): # TODO refactor this function
        for p in self.emof_properties:
            f = self.emof_functions.get(p["measurementType"])
            if f:
                mv, mvid = f(
                    observatory=self.observatory,
                    omics=self.omics,
                    sampling=self.sampling,
                    pema=self.pema,
                    ptax=self.ptax if self.genomic_region == "COI" else None
                )
            else:
                mv = p["measurementValue"]
                mvid = p["measurementValueID"]
            self.extended_measurement_or_facts_extension.add_node(
                occurrenceID=self.occurrence_id,  # node identifier
                measurementType=p["measurementType"],
                measurementUnit=p["measurementUnit"],
                measurementValue=mv,
                measurementTypeID=p["measurementTypeID"],
                measurementUnitID=p["measurementUnitID"],
                measurementValueID=mvid,
            )

    def _serialize(self):
        self.occurrence_core.serialize()
        self.dna_extension.serialize()
        self.extended_measurement_or_facts_extension.serialize()
