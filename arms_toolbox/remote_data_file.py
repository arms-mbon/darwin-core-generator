import requests
from Bio import SeqIO
from io import StringIO
import pandas as pd


class RemoteDataFile:
    def __init__(self):
        self.download_url = None

    def fetch(self):
        return NotImplementedError


class RemoteExcelFile(RemoteDataFile):
    def fetch(self, csv_fallback=False):
        if csv_fallback:
            try:
                return pd.read_excel(self.download_url)
            except:
                return pd.read_csv(self.download_url.replace(".xlsx", ".csv"))
        else:
            return pd.read_excel(self.download_url)


class RemoteCSVFile(RemoteDataFile):
    def fetch(self):
        return pd.read_csv(self.download_url)


class RemoteTSVFile(RemoteDataFile):
    def fetch(self):
        ...


class ObservatoryDataFrame(RemoteCSVFile):
    def __init__(self):
        self.download_url = "https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_ObservatoryData.csv"
    
        
class OmicsDataFrame(RemoteCSVFile):
    def __init__(self):
        self.download_url = "https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_OmicsData.csv"


class SamplingDataFrame(RemoteCSVFile):
    def __init__(self):
        self.download_url = "https://raw.githubusercontent.com/arms-mbon/data_workspace/main/qualitycontrolled_data/combined/combined_SamplingEventData.csv"


class PEMADataFrame(RemoteExcelFile):
    def __init__(self, genomic_region, time_window, data_release):
        self.genomic_region = genomic_region
        self.time_window = time_window
        self.data_release = data_release

        if self.genomic_region == "ITS":
            self.download_url = f"https://raw.githubusercontent.com/arms-mbon/analysis_release_{self.data_release}/main/taxonomic_assignments/Extended_final_table_{self.time_window}_{self.genomic_region}_noBlank.xlsx"
        elif self.genomic_region == "18S":
            self.download_url = f"https://raw.githubusercontent.com/arms-mbon/analysis_release_{self.data_release}/main/taxonomic_assignments/Extended_final_table_{self.time_window}_{self.genomic_region}_noBlank_TaxonomyCurated.xlsx"
        elif self.genomic_region == "COI":
            self.download_url = f"https://raw.githubusercontent.com/arms-mbon/analysis_release_{self.data_release}/main/taxonomic_assignments/Extended_final_table_{self.time_window}_{self.genomic_region}_noBlank_TaxonomyFull.csv"
        else:
            raise AssertionError
    
    def fetch(self):
        df = super().fetch(csv_fallback=True)
        if self.genomic_region == "COI":
            return df.rename(columns={"ASV_number:amplicon": "OTU"})            
        else:
            return df


class PTAXDataFrame(RemoteTSVFile):
    """PEMA Taxonomy
    """
    def __init__(self, time_window, data_release):
        self.time_window = time_window
        self.data_release = data_release
        self.download_url = f"https://raw.githubusercontent.com/arms-mbon/analysis_release_{self.data_release}/main/taxonomic_assignments/tax_assignments_{self.time_window}_COI_noBlank.tsv"

    def fetch(self):
        df = pd.read_csv(self.download_url, sep=r"\s+", header=None)
        df[0] = df[0].apply(lambda _: _.split("_")[0])
        df.columns = [
            "amplicon",
            "kingdom", "kingdom_confidence",
            "phylum", "phylum_confidence",
            "class", "class_confidence",
            "order", "order_confidence",
            "g",
            "family", "family_confidence",
            "genus", "genus_confidence",
            "species", "species_confidence",
        ]
        return df


class MDAFasta(RemoteDataFile):
    def __init__(self, genomic_region, download_url):
        self.genomic_region = genomic_region
        self.download_url = download_url

    def fetch(self):
        response = requests.get(self.download_url.replace("directlink.php?fid", "download.php?file")).content.decode()
        seq = SeqIO.to_dict(SeqIO.parse(StringIO(response), "fasta"))
        if self.genomic_region == "COI":  # grouped fasta
            seq_split = {k.split("_")[0]: v for k, v in seq.items()}
            assert len(seq_split) == len(seq)
            return seq_split
        else:
            return seq
