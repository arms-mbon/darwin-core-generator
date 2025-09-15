import contextlib
import requests
import time
import pandas as pd
from io import StringIO
from functools import cache

@cache
def get_aphia_taxonomy(ncbi_tax_id, sleep=1):
    """Get the Aphia ID for the NCBI ID"""
    response = requests.get(
        f"https://www.marinespecies.org/rest/AphiaRecordByExternalID/{ncbi_tax_id}?type=ncbi",
        headers={"Accept": "application/json"}
    )
    time.sleep(sleep)
    if response.status_code == 200:
        return response.json()
    else:
        return {}

@cache
def get_col_scientific_name(ncbi_scientific_name, sleep=1):
    response = requests.get(
        f"https://api.checklistbank.org/nidx/match?q={ncbi_scientific_name}",
        headers={"Accept": "application/json"}
    )
    time.sleep(sleep)

    if response.status_code == 200:
        body = response.json()
    else:
        body = {}
    
    if "name" in body:
        return body["name"].get("scientificName", "")
    else:
        return ""

@cache
def get_ena_taxonomy(ncbi_tax_id, sleep=1):
    """ Get the rank information for the species name from ENA"""
    response = requests.get(
        f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{ncbi_tax_id}?binomialOnly=false",
        headers={"Accept": "application/json"}
    )
    time.sleep(sleep)
    if response.status_code == 200:
        return response.json()
    else:
        return {}

@cache
def get_ena_sample_accession(run_accessions, sleep=1):
    """Get the sample accession number for a specific set of run accession numbers"""
    response = requests.post(
        "https://www.ebi.ac.uk/ena/portal/api/search",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data={
            "result": "read_run",
            "fields": "sample_accession",
            "includeAccessionType": "run",
            "includeAccessions": run_accessions,
            "format": "tsv",
        }
    )
    time.sleep(sleep)
    if response.status_code == 200:
        return pd.read_csv(StringIO(response.text), sep="\t")["sample_accession"].iloc[0]
    else:
        return "None"

@cache
def get_ncbi_tax_id(scientific_name, sleep=1):
    """Get the NCBI ID for a species name"""
    with contextlib.redirect_stdout(StringIO()):
        from arms_toolbox.ncbi_taxonomist_wrapper import resolve
        from copy import deepcopy
        import sys
        sys_exit_copy = deepcopy(sys.exit)
        sys.exit = lambda *args, **kwargs: None
        ncbi_tax_id = resolve(scientific_name)
        sys.exit = sys_exit_copy
        time.sleep(sleep)
        return ncbi_tax_id
