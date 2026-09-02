# `arms_toolbox`

> **Disclaimer regarding the use of Generative AI in this document**
>
> This documentation was generated with the assistance of a large language model (LLM) based on inspection of the `arms_toolbox` source code and repository documentation. It has been reviewed by a developer for technical accuracy and consistency with the current implementation. The source code remains the authoritative reference for package behavior.



`arms_toolbox` is the Python package underlying the ARMS-MBON `darwin-core-generator`. It transforms ARMS-MBON sampling metadata, omics metadata, PEMA taxonomic assignments, and amplicon FASTA files into Darwin Core-compatible CSV datasets. The current implementation supports the **18S**, **ITS**, and **COI** genomic regions (i.e. marker genes) and generates an Occurrence Core, DNA-derived data extension, and Extended Measurement or Facts (eMoF) extension.

## Architecture

The package follows a lightweight ETL architecture:

```text
ARMS-MBON repositories / MDA
            │
            ▼
         extract
            │
            ▼
        transform
       ┌────┼────┐
       │    │    │
      18S  ITS  COI
       │    │    │
       └────┼────┘
            ▼
           load
            │
            ▼
      occurrence.csv
      dnaextension.csv
      emof.csv
```

`Pipeline` provides the common workflow. `Pipeline18S`, `PipelineITS`, and `PipelineCOI` configure this base class for the supported marker regions. The COI implementation additionally enables COI-specific eMoF properties and transformation functions.

## Installation
Runtime dependencies are maintained in the repository-level `requirements.txt` and currently include:

```text
pandas
biopython
pycountry
ncbi-taxonomist
requests
openpyxl
tqdm
pyyaml
```

These dependencies can be installed using `pip`:
```
pip install -r requirements.txt
```

`arms_toolbox` currently has no independent packaging metadata such as `pyproject.toml` or `setup.py`. It should therefore be treated as a source package inside the repository rather than something installed independently with `pip install arms_toolbox`.

## Running a pipeline

A minimal invocation is, for example:

```python
from arms_toolbox.pipeline_18s import Pipeline18S

pipeline = Pipeline18S(
    time_window="April2021",
    aligned_assignment_url="<MDA FASTA URL>",
    dataset_name="data_release_002",
    data_release="002",
)

pipeline.run()
```

## Inputs

A pipeline run depends on several separately maintained datasets:

```text
data_workspace
├── combined_ObservatoryData.csv
├── combined_OmicsData.csv
└── combined_SamplingEventData.csv

analysis_release_<N>
├── taxonomic assignment table
└── COI PEMA tax_assignments TSV

VLIZ Marine Data Archive
└── aligned/grouped FASTA file

darwin-core-generator/data/schemas
└── Darwin Core schema files
```

## Outputs

For each `<time_window>` and `<genomic_region>` pair, the package generates:

```text
data/output/<time_window>_<genomic_region>/occurrence.csv
data/output/<time_window>_<genomic_region>/dnaextension.csv
data/output/<time_window>_<genomic_region>/emof.csv
```

## Error handling

Extraction and transformation exceptions are logged and processing continues, rather than aborting the entire pipeline run. This behavior is suitable for batch processing but means successful pipeline termination does not necessarily imply that every input occurrence was added to the output datasets. Consumers should inspect application logs and validate record counts after execution.

## Package layout in-depth

```text
arms_toolbox/
├── __init__.py
├── data_sink.py
├── data_source.py
├── function.py
├── graphoid.py
├── ncbi_taxonomist_wrapper.py
├── parameter.py
├── pipeline.py
├── pipeline_18s.py
├── pipeline_coi.py
└── pipeline_its.py
```

### `pipeline.py`

Defines the main `Pipeline` abstraction:

```python
Pipeline(
    genomic_region,
    time_window,
    aligned_assignment_url,
    dataset_name,
    data_release,
)
```

The constructor records the marker region, sampling time window, MDA FASTA URL, target dataset, and ARMS-MBON data-release identifier. It also installs the mapping between eMoF properties and the transformation functions defined in `function.py`. During execution the pipeline loads its remote data sources, iterates through replicate material sample IDs, joins the relevant sampling and omics records, iterates over taxonomic assignments, creates an occurrence identifier, transforms the records into Darwin Core fields, and finally serializes the output tables. Extraction or transformation failures are logged and processing continues with subsequent records.

The occurrence identifier is constructed from the replicate material sample identifier and the PEMA OTU identifier:

```python
occurrence_id = f"{ReplicateMaterialSampleID}:{OTU}"
```

Taxonomic records where both the classification and NCBI taxon identifier are reported as `not found` are rejected. Sampling metadata are also required to contain exactly one record for each replicate material sample ID being processed.

### Marker-specific pipelines

`Pipeline18S` fixes `genomic_region="18S"`. `PipelineITS` behaves equivalently with `genomic_region="ITS"`. `PipelineCOI` fixes `genomic_region="COI"` and extends the eMoF transformation map with `originalScientificName`, `originalTaxonRank`, and `originalScientificNameConfidenceLevel`.

### `data_source.py`

Contains adapters for retrieving the pipeline's input datasets. The generic classes are `RemoteDataFile`, `RemoteExcelFile`, `RemoteCSVFile`, and the currently unimplemented `RemoteTSVFile`. Concrete data sources include `ObservatoryDataFrame`, `OmicsDataFrame`, `SamplingDataFrame`, `PEMADataFrame`, `PTAXDataFrame`, and `MDAFasta`.

The observatory, omics, and sampling-event datasets are read from the combined quality-controlled files in the ARMS-MBON `data_workspace` repository. PEMA assignments are read from the corresponding `analysis_release_<release>` repository. File naming varies by marker. ITS reads `Extended_final_table_<time_window>_ITS_noBlank.xlsx`. 18S reads `Extended_final_table_<time_window>_18S_noBlank_TaxonomyCurated.xlsx`. COI reads `Extended_final_table_<time_window>_COI_noBlank_TaxonomyFull.csv` and renames `ASV_number:amplicon` to `OTU`.

COI also loads a `tax_assignments_<time_window>_COI_noBlank.tsv` PEMA taxonomy file. This is parsed into taxonomic ranks and associated confidence values.

`MDAFasta` retrieves sequence data from a VLIZ Marine Data Archive download link and parses FASTA records with Biopython. COI uses grouped FASTA identifiers and normalizes each identifier by taking the component preceding the first underscore.

### `data_sink.py`

Defines three schema-backed outputs:

```python
OccurrenceCore(genomic_region, time_window)
DNAExtension(genomic_region, time_window)
ExtendedMeasurementOrFactsExtension(genomic_region, time_window)
```

The corresponding files are written beneath:

```text
data/output/<time_window>_<genomic_region>/
├── occurrence.csv
├── dnaextension.csv
└── emof.csv
```

Occurrence and DNA-extension schemas are marker-specific. The eMoF table uses the shared `extended_measurement_or_facts_extension_schema.json`.

### `graphoid.py`

`Graphoid` is the internal schema-backed table abstraction used by the output classes. Despite its graph-oriented interface, its current backing implementation is a Pandas `DataFrame`. A JSON schema determines the available columns and their default values. `add_node()` appends a row and `serialize()` creates the destination directory before writing the DataFrame as CSV. The code describes this abstraction as an intermediate step toward a possible object-graph mapper (cf. [py-sema's OGM](https://github.com/vliz-be-opsci/py-sema/tree/main/sema/commons/ogm) & [Darwin Core RDF Guide](https://dwc.tdwg.org/rdf/)).

```python
graph = Graphoid("output.csv", "schema.json")
graph.add_node(occurrenceID="sample:ASV1")
graph.serialize()
```

### `function.py`

Contains the transformation functions used to populate Darwin Core and eMoF values. `Pipeline` registers these functions against property names and invokes them while constructing output rows. Consequently, most domain-specific conversion logic belongs here rather than in the marker-specific pipeline classes.

When adding a new derived field, the recommended procedure is to implement its transformation in `function.py`, add the appropriate schema/property definition under `data/schemas`, and register the function in `Pipeline.emof_functions` or the appropriate marker subclass.

### `parameter.py`

Loads eMoF property definitions when the module is imported:

```python
EMOFProperties
EMOFPropertiesCOI
```

These originate from JSON files under `data/schemas`.

### `ncbi_taxonomist_wrapper.py`

Provides:

```python
resolve(name)
```

This function wraps `ncbi-taxonomist` to remotely resolve a taxonomic name to an NCBI taxonomy identifier. Resolution errors are converted to `None`. Wrapping is required to circumvent the use of the `ncbi-taxonomist` CLI.
