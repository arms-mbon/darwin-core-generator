from .pipeline_its import PipelineITS

class Pipeline18S(PipelineITS):
    def __init__(self, time_window, aligned_assignment_url):
        super().__init__(time_window, aligned_assignment_url)
        self.genomic_region = "18S"
        self.occurrence_core_schema = "./data/schemas/occurrence_core_schema_18s.json"
        self.dna_extension_schema = "./data/schemas/dna_extension_schema_18s.json"
