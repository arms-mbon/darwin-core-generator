from .pipeline import Pipeline

class PipelineITS(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(genomic_region="ITS", *args, **kwargs)
