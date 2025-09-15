from datetime import datetime
from pycountry import countries
from .remote_data_api import get_ena_sample_accession, get_col_scientific_name, get_aphia_taxonomy, get_ena_taxonomy, get_ncbi_tax_id
from functools import cache

def _parse_node_label(taxon_ncbi_tax_id):
    # print(taxon_ncbi_tax_id)
    taxon = taxon_ncbi_tax_id.split(":")[0]
    ncbi_tax_ids = taxon_ncbi_tax_id.split(":")[1].split("\\n")
    ncbi_tax_id = ncbi_tax_ids[0] if ncbi_tax_ids else None
    return taxon, ncbi_tax_id

def _preprocess_classification(classification):
    prefix = "Cellular organisms;"
    if classification.startswith(prefix):
        classification = classification[len(prefix):]
    return classification.replace("_", " ")

def _parse_classification(classification, taxon):
    lineage = (
        _preprocess_classification(classification)
        .split(taxon)[0][:-1]
        .split(";")
    )
    return lineage

@cache
def _walk_phylogeny(node_label=None, taxon=None, ncbi_tax_id=None, classification=None, lineage=None):
    assert node_label or (taxon or ncbi_tax_id)
    if node_label:
        taxon, ncbi_tax_id = _parse_node_label(node_label)
    
    if not ncbi_tax_id:
        assert taxon
        ncbi_tax_id = get_ncbi_tax_id(taxon)

    if ncbi_tax_id:
        atax = get_aphia_taxonomy(ncbi_tax_id)
        if atax:
            aphia_id = atax["AphiaID"]
            scientific_name_id = f"urn:lsid:marinespecies.org:taxname:{aphia_id}"
            scientific_name = atax["scientificname"]
            taxon_rank = atax["rank"]
            return scientific_name_id, scientific_name, taxon_rank

        etax = get_ena_taxonomy(ncbi_tax_id)
        if not taxon:
            taxon = etax["scientificName"]
        col_scientific_name = get_col_scientific_name(taxon)
        if taxon == col_scientific_name:
            scientific_name_id = ""
            scientific_name = taxon
            taxon_rank = etax["rank"]
            return scientific_name_id, scientific_name, taxon_rank
    
    if classification:
        lineage = _parse_classification(classification, taxon)

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

def submerged_time(**kwargs):
    date_deployed = kwargs["df_sampling_subset"]["DateDeployed"].iloc[0]
    date_collected = kwargs["df_sampling_subset"]["DateCollected"].iloc[0]
    mv = (
        datetime.strptime(date_collected, "%Y-%m-%d")
        - datetime.strptime(date_deployed, "%Y-%m-%d")
    ).days
    return mv, "" 

def preservative_used(**kwargs):
    mv = kwargs["df_sampling_subset"]["Preservative"].iloc[0]
    mvid = {
        "etoh": "https://pubchem.ncbi.nlm.nih.gov/compound/702",
        "dmso": "https://pubchem.ncbi.nlm.nih.gov/compound/679"
    }.get(mv.lower(), "")
    return mv, mvid

def lower_limit_filter_size(**kwargs):
    mv = kwargs["df_sampling_subset"]["Filter"].iloc[0]
    return mv, "" 

def field_replicate(**kwargs):
    mv = kwargs["df_sampling_subset"]["FieldReplicate"].iloc[0]
    return mv, "" 

def technical_replicate(**kwargs):
    if (kwargs["df_sampling_subset"]["SequencingRunRepeat"].iloc[0] == "second sequencing run (repeat)"):
        mv = kwargs["df_sampling_subset"]["MaterialSampleID"].iloc[0]
        return mv, "" 
    else:
        return "", "" 

def original_scientific_name(**kwargs):
    row = kwargs["row"]
    ptax = kwargs["ptax"]
    otu_id = row["OTU"].split(":")[1]
    try:
        matchrow = ptax.loc[ptax["amplicon"] == otu_id]
        species_name = matchrow["species"].values[0]
        _ = matchrow["species_confidence"].values[0]
        return species_name, ""
    except Exception:
        return "none", ""

def original_taxon_rank(**kwargs):
    row = kwargs["row"]
    ptax = kwargs["ptax"]
    otu_id = row["OTU"].split(":")[1]
    try:
        matchrow = ptax.loc[ptax["amplicon"] == otu_id]
        _ = matchrow["species"].values[0]
        _ = matchrow["species_confidence"].values[0]
        return "species", ""
    except Exception:
        return "none", ""

def original_scientific_name_confidence_level(**kwargs):
    row = kwargs["row"]
    ptax = kwargs["ptax"]
    otu_id = row["OTU"].split(":")[1]
    try:
        matchrow = ptax.loc[ptax["amplicon"] == otu_id]
        _ = matchrow["species"].values[0]
        conflevel = matchrow["species_confidence"].values[0]
        return conflevel, ""
    except Exception:
        return "0", ""



def material_sample_id(associated_sequences):
    return get_ena_sample_accession(associated_sequences)

def country_code(df):
    return countries.get(alpha_3=df["Country"].iloc[0]).alpha_2

def habitat(df):
    return "; ".join([
        df["Monitoring area"].iloc[0].strip(),
        df["Anthropogenic influence"].iloc[0].strip(),
        df["IUCN habitat type"].iloc[0].strip(),
    ])

def occurrence_remarks(node_label):
    _, ncbi_tax_id = _parse_node_label(node_label)
    atax = get_aphia_taxonomy(ncbi_tax_id)
    if str(atax.get("isTerrestrial")) == "1":
        return "uncertain classification (not marine)"
    else:
        return ""

def event_remarks(df):
    crate_cover = df["CrateCover"].iloc[0]
    fraction = df["Fraction"].iloc[0]
    return f"Was a crate cover used when retrieving the unit: {crate_cover}; Which component of the ARMS unit material was extracted: {fraction}"

def verbatim_identification(classification):
    return _preprocess_classification(classification)

def scientific_name_id(node_label, classification):
    scientific_name_id, _, _ = _walk_phylogeny(
        node_label=node_label,
        classification=classification,
    )
    return scientific_name_id

def scientific_name(node_label, classification):
    _, scientific_name, _ = _walk_phylogeny(
        node_label=node_label,
        classification=classification,
    )
    return scientific_name 

def taxon_rank(node_label, classification):
    _, _, taxon_rank = _walk_phylogeny(
        node_label=node_label,
        classification=classification,
    )
    return taxon_rank

def dna_sequence(otu, fasta, grouped=False):
    if grouped:
        return fasta[otu.split(":")[1]].seq
    else:
        return fasta[otu].seq

