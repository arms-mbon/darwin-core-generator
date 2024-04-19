from pathlib import Path

class DwCArchive:
    def __init__(self):
        self.core = None
        self.extensions = []
    
    def read(self):
        ...

    def write(self, path, compression=None):
        Path(path).mkdir(parents=True, exist_ok=True)
        tables = [self.core] + self.extensions
        for table in tables:
            table.to_csv(root_path=path)
    
    def add_core(self, core, file_name=None):
        self.assign_file_name(core, file_name)
        self.core = core
    
    def add_extension(self, extension, file_name=None):
        self.assign_file_name(extension, file_name)
        self.extensions.append(extension)
    
    def add_eml(self, eml):
        """
        eml is either a path to an eml file or an xml string
        """
        ...

    @staticmethod
    def assign_file_name(table, file_name):
        if not file_name:
            table.file_name = f"{type(table).__name__.lower()}.txt"
        else:
            table.file_name = file_name
