from .pipeline import Pipeline
from .remote_data_file import PTAXDataFrame
import arms_toolbox.function as fn
import json

class PipelineCOI(Pipeline):
    def __init__(self, time_window, data_release, *args, **kwargs):
        super().__init__(
            genomic_region="COI",
            time_window=time_window,
            data_release=data_release,
            *args,
            **kwargs
        )

        # input
        self.ptax = PTAXDataFrame(time_window, data_release).fetch()

        # output
        self.emof_properties = json.load(open("./data/schemas/extended_measurement_or_facts_extension_properties_schema_coi.json"))
        self.emof_functions.update({
            "originalScientificName": fn.original_scientific_name,
            "originalTaxonRank": fn.original_taxon_rank,
            "originalScientificNameConfidenceLevel": fn.original_scientific_name_confidence_level,
        })

