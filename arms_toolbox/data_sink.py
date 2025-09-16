from .graphoid import Graphoid


class OccurrenceCore(Graphoid):
    def __init__(self, genomic_region, time_window):
        self.sink_path=f"./data/output/{time_window}_{genomic_region}/occurrence.csv"
        self.schema_path=f"./data/schemas/occurrence_core_schema_{genomic_region.lower()}.json"
        super().__init__(self.sink_path, self.schema_path)


class DNAExtension(Graphoid):
    def __init__(self, genomic_region, time_window):
        self.sink_path=f"./data/output/{time_window}_{genomic_region}/dnaextension.csv"
        self.schema_path=f"./data/schemas/dna_extension_schema_{genomic_region.lower()}.json"
        super().__init__(self.sink_path, self.schema_path)


class ExtendedMeasurementOrFactsExtension(Graphoid):
    def __init__(self, genomic_region, time_window):
        self.sink_path=f"./data/output/{time_window}_{genomic_region}/emof.csv"
        self.schema_path="./data/schemas/extended_measurement_or_facts_extension_schema.json"
        super().__init__(self.sink_path, self.schema_path)
