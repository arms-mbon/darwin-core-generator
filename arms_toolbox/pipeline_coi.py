from .pipeline import Pipeline
from .parameter import EMOFPropertiesCOI
import arms_toolbox.function as fn

class PipelineCOI(Pipeline):
    def __init__(self, time_window, data_release, *args, **kwargs):
        super().__init__(
            genomic_region="COI",
            time_window=time_window,
            data_release=data_release,
            *args,
            **kwargs
        )

        # parameters
        self.emof_properties = EMOFPropertiesCOI  # TODO refactor this
        self.emof_functions.update({  # TODO refactor this
            "originalScientificName": fn.original_scientific_name,
            "originalTaxonRank": fn.original_taxon_rank,
            "originalScientificNameConfidenceLevel": fn.original_scientific_name_confidence_level,
        })
