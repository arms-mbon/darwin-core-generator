from .pipeline import Pipeline

class PipelineCOI(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(genomic_region="COI", *args, **kwargs)
