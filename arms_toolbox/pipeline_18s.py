from .pipeline import Pipeline

class Pipeline18S(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(genomic_region="18S", *args, **kwargs)
