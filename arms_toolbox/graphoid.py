"""
TODO: end goal is to work with an object-graph mapper. Graphoid is already an intermediate, as it provides an ogm-like interface, but still is a dataframe internally.
"""
import json
import pandas as pd
from pathlib import Path

class Graphoid:
    def __init__(self, sink_path, schema_path):
        self.sink_path = Path(sink_path)
        self.schema = json.load(open(schema_path))  # TODO replace json by yaml
        self.df = pd.DataFrame(columns=self.schema.keys())
    
    def add_node(self, **kwargs):
        row = self.schema.copy()
        for k, v in kwargs.items():
            row[k] = v
        self.df.loc[len(self.df)] = row.values()

    def serialize(self):
        self.sink_path.parent.mkdir(parents=True, exist_ok=True)
        self.df.to_csv(self.sink_path, index=False)
