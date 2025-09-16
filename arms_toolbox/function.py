import contextlib
import requests
import time
import pandas as pd
import logging
from io import StringIO
from functools import cache
from datetime import datetime
from pycountry import countries
from functools import cache

logger = logging.getLogger(__name__)

def _parse_pema_leaf(pema_leaf):
    taxon = pema_leaf.split(":")[0]
    ncbi_tax_ids = pema_leaf.split(":")[1].split("\\n")
    ncbi_tax_id = ncbi_tax_ids[0] if ncbi_tax_ids else None
    return taxon, ncbi_tax_id

def _preprocess_pema_tree(pema_tree):
    prefix = "Cellular organisms;"
    if pema_tree.startswith(prefix):
        pema_tree = pema_tree[len(prefix):]
    return pema_tree.replace("_", " ")

def _parse_pema_tree(classification, taxon):
    lineage = (
        _preprocess_pema_tree(classification)
        .split(taxon)[0][:-1]
        .split(";")
    )
    return lineage

@cache
def _walk_phylogeny(pema_leaf=None, taxon=None, ncbi_tax_id=None, pema_tree=None, lineage=None):
    assert pema_leaf or (taxon or ncbi_tax_id)
    if pema_leaf:
        taxon, ncbi_tax_id = _parse_pema_leaf(pema_leaf)
    
    if not ncbi_tax_id:
        assert taxon
        ncbi_tax_id = _get_ncbi_tax_id(taxon)

    if ncbi_tax_id:
        atax = _get_aphia_taxonomy(ncbi_tax_id)
        if atax:
            aphia_id = atax["AphiaID"]
            scientific_name_id = f"urn:lsid:marinespecies.org:taxname:{aphia_id}"
            scientific_name = atax["scientificname"]
            taxon_rank = atax["rank"]
            return scientific_name_id, scientific_name, taxon_rank

        etax = _get_ena_taxonomy(ncbi_tax_id)
        if not taxon:
            taxon = etax["scientificName"]
        col_scientific_name = _get_col_scientific_name(taxon)
        if taxon == col_scientific_name:
            scientific_name_id = ""
            scientific_name = taxon
            taxon_rank = etax["rank"]
            return scientific_name_id, scientific_name, taxon_rank
    
    if pema_tree:
        assert taxon
        lineage = _parse_pema_tree(pema_tree, taxon)

    if lineage:
        if isinstance(lineage, str):
            lineage_unwrapped = lineage.split(";")
        else:
            lineage_unwrapped = lineage
        taxon = lineage_unwrapped[-1]
        lineage_wrapped = ";".join(lineage_unwrapped[:-1])
        return _walk_phylogeny(taxon=taxon, lineage=lineage_wrapped)
    else:
        scientific_name_id = "urn:lsid:marinespecies.org:taxname:1"
        scientific_name = "Biota"
        taxon_rank = ""
        return scientific_name_id, scientific_name, taxon_rank

@cache
def _get_aphia_taxonomy(ncbi_tax_id, sleep=1):
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
def _get_col_scientific_name(ncbi_scientific_name, sleep=1):
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
def _get_ena_taxonomy(ncbi_tax_id, sleep=1):
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
def _get_ena_sample_accession(run_accessions, sleep=1):
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
def _get_ncbi_tax_id(scientific_name, sleep=1):
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

def submerged_time(**kwargs):
    date_deployed = kwargs["sampling"]["DateDeployed"]
    date_collected = kwargs["sampling"]["DateCollected"]
    mv = (
        datetime.strptime(date_collected, "%Y-%m-%d")
        - datetime.strptime(date_deployed, "%Y-%m-%d")
    ).days
    return mv, "" 

def preservative_used(**kwargs):
    mv = kwargs["sampling"]["Preservative"]
    mvid = {
        "etoh": "https://pubchem.ncbi.nlm.nih.gov/compound/702",
        "dmso": "https://pubchem.ncbi.nlm.nih.gov/compound/679"
    }.get(mv.lower(), "")
    return mv, mvid

def lower_limit_filter_size(**kwargs):
    mv = kwargs["sampling"]["Filter"]
    return mv, "" 

def field_replicate(**kwargs):
    mv = kwargs["sampling"]["FieldReplicate"]
    return mv, "" 

def technical_replicate(**kwargs):
    if (kwargs["sampling"]["SequencingRunRepeat"] == "second sequencing run (repeat)"):
        mv = kwargs["sampling"]["MaterialSampleID"]
        return mv, "" 
    else:
        return "", "" 

def original_scientific_name(**kwargs):
    pema = kwargs["pema"]
    ptax = kwargs["ptax"]
    otu_id = pema["OTU"].split(":")[1]
    try:
        matchrow = ptax.loc[ptax["amplicon"] == otu_id]
        species_name = matchrow["species"].values[0]
        _ = matchrow["species_confidence"].values[0]
        return species_name, ""
    except Exception:
        return "none", ""

def original_taxon_rank(**kwargs):
    pema = kwargs["pema"]
    ptax = kwargs["ptax"]
    otu_id = pema["OTU"].split(":")[1]
    try:
        matchrow = ptax.loc[ptax["amplicon"] == otu_id]
        _ = matchrow["species"].values[0]
        _ = matchrow["species_confidence"].values[0]
        return "species", ""
    except Exception:
        return "none", ""

def original_scientific_name_confidence_level(**kwargs):
    pema = kwargs["pema"]
    ptax = kwargs["ptax"]
    otu_id = pema["OTU"].split(":")[1]
    try:
        matchrow = ptax.loc[ptax["amplicon"] == otu_id]
        _ = matchrow["species"].values[0]
        conflevel = matchrow["species_confidence"].values[0]
        return conflevel, ""
    except Exception:
        return "0", ""

def material_sample_id(associated_sequences):
    return _get_ena_sample_accession(associated_sequences)

def country_code(observatory: pd.Series):
    return countries.get(alpha_3=observatory["Country"]).alpha_2

def habitat(observatory: pd.Series):
    return "; ".join([
        observatory["Monitoring area"].strip(),
        observatory["Anthropogenic influence"].strip(),
        observatory["IUCN habitat type"].strip(),
    ])

def occurrence_remarks(pema_leaf):
    _, ncbi_tax_id = _parse_pema_leaf(pema_leaf)
    atax = _get_aphia_taxonomy(ncbi_tax_id)
    if str(atax.get("isTerrestrial")) == "1":
        return "uncertain classification (not marine)"
    else:
        return ""

def event_remarks(sampling: pd.Series):
    crate_cover = sampling["CrateCover"]
    fraction = sampling["Fraction"]
    return f"Was a crate cover used when retrieving the unit: {crate_cover}; Which component of the ARMS unit material was extracted: {fraction}"

def verbatim_identification(pema_tree):
    return _preprocess_pema_tree(pema_tree)

def scientific_name_id(pema_leaf, pema_tree):
    scientific_name_id, _, _ = _walk_phylogeny(
        pema_leaf=pema_leaf,
        pema_tree=pema_tree,
    )
    return scientific_name_id

def scientific_name(pema_leaf, pema_tree):
    _, scientific_name, _ = _walk_phylogeny(
        pema_leaf=pema_leaf,
        pema_tree=pema_tree,
    )
    return scientific_name 

def taxon_rank(pema_leaf, pema_tree):
    _, _, taxon_rank = _walk_phylogeny(
        pema_leaf=pema_leaf,
        pema_tree=pema_tree,
    )
    return taxon_rank

def dna_sequence(otu, mda_fasta, genomic_region):
    if genomic_region == "COI": # grouped
        return mda_fasta[otu.split(":")[1]].seq
    else:
        return mda_fasta[otu].seq

