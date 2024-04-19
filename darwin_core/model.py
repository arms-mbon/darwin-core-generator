import json
import pandas as pd
from pathlib import Path

class Table:
    def __init__(self, schema_path):
        self.schema = json.load(open(schema_path))
        self.df = pd.DataFrame(columns=self.schema.keys())
        self.file_name = None
    
    def add_row(self, **kwargs):
        row = self.schema.copy()
        for k, v in kwargs.items():
            row[k] = v
        self.df.loc[len(self.df)] = row.values()

    def to_csv(self, root_path=None):
        if root_path:
            file_path = Path(root_path) / self.file_name
        else:
            file_path = self.file_name
        self.df.to_csv(file_path, index=False)

class Core(Table):
    ...

class Extension(Table):
    ...

class OccurrenceCore(Core):
    ...

class DNAExtension(Extension):
    ...

class ExtendedMeasurementOrFactsExtension(Extension):
    ...
