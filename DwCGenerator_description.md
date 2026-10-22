# ARMS-MBON Darwin Core Generator — Pipeline Overview \& File Reference

## 1\. Project Overview

The **Darwin Core Generator** (`darwin-core-generator`) is an Extract-Transform-Load (ETL) toolkit built for the **ARMS-MBON** (Autonomous Reef Monitoring Structures – Marine Biodiversity Observation Network) project. Its primary goal is to transform raw metabarcoding bioinformatics outputs from the **PEMA** pipeline—across three genetic markers (**18S**, **ITS**, and **COI**)—along with field sampling, observatory, and omics metadata into standardized **Darwin Core Archive (DwC-A)** CSV tables ready for submission to **EurOBIS** / **OBIS**.

For each combination of **sampling time window** (e.g., `September2020`, `April2021`, `August2023`) and **genomic region** (`18S`, `ITS`, `COI`), the pipeline produces three linked tables joined by `occurrenceID`:

1. **`occurrence.csv`** — Occurrence Core record (who, what, where, when, how many reads).
2. **`dnaextension.csv`** — DNA-derived Data Extension (amplicon sequence, primers, SOPs, ENVO habitat terms).
3. **`emof.csv`** — Extended Measurement or Facts (eMoF) Extension (submerged time, preservative, filter size, replicates, and COI original taxonomic assignment metrics).

\---

## 2\. End-to-End Pipeline Architecture

The core workflow is orchestrated by [`Pipeline`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline.py#L14-L177) in [pipeline.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline.py) and specialized by marker-specific subclasses ([`Pipeline18S`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_18s.py#L3-L5), [`PipelineITS`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_its.py#L3-L5), and [`PipelineCOI`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_coi.py#L5-L21)).

```mermaid
flowchart TB
    subgraph Sources\["1. Remote Data Sources (Extract)"]
        DW\["GitHub: data\_workspace<br/>• combined\_ObservatoryData.csv<br/>• combined\_OmicsData.csv<br/>• combined\_SamplingEventData.csv"]
        AR\["GitHub: analysis\_release\_002<br/>• Extended\_final\_table (18S/ITS/COI)<br/>• tax\_assignments TSV (COI only)"]
        MDA\["VLIZ Marine Data Archive (MDA)<br/>• Aligned / Grouped FASTA files"]
        SCH\["Local JSON Schemas<br/>• data/schemas/\*.json"]
    end

    subgraph Engine\["2. Core ETL Engine (arms\_toolbox)"]
        INIT\["Initialize Pipeline<br/>Fetch remote tables \& FASTA into memory"]
        LOOP\_RMSI\["Iterate over ReplicateMaterialSampleID (rmsi) from extended\_final\_table (omics)<br/>Validate unique match in Sampling, Omics \& Observatory"]
        FILTER\["Filter PEMA OTUs / ASVs<br/>Keep rows where read abundance > 1<br/>Reject 'not found' classifications"]
        OCC\_ID\["Construct Key<br/>occurrenceID = '{rmsi}:{OTU}'"]
    end

    subgraph APIs\["3. External Taxonomic \& Accession APIs (Transform)"]
        ENA\_ACC\["ENA Portal API<br/>Run Accession → Sample Accession"]
        WORMS\["WoRMS REST API<br/>NCBI TaxID → AphiaID \& Marine Check"]
        ENA\_TAX\["ENA Taxonomy API<br/>NCBI TaxID → Rank \& Scientific Name"]
        COL\["Catalogue of Life (ChecklistBank)<br/>Verify Official Scientific Name"]
        NCBITAX\["ncbi-taxonomist<br/>Taxon Name → NCBI TaxID"]
    end

    subgraph Sinks\["4. Schema-Backed Sinks (Load / Serialize)"]
        G\_OCC\["OccurrenceCore (Graphoid)<br/>→ occurrence.csv"]
        G\_DNA\["DNAExtension (Graphoid)<br/>→ dnaextension.csv"]
        G\_EMOF\["ExtendedMeasurementOrFacts (Graphoid)<br/>→ emof.csv"]
    end

    DW \& AR \& MDA \& SCH --> INIT
    INIT --> LOOP\_RMSI --> FILTER --> OCC\_ID
    OCC\_ID <--> ENA\_ACC \& WORMS \& ENA\_TAX \& COL \& NCBITAX
    OCC\_ID --> G\_OCC \& G\_DNA \& G\_EMOF
```

### Step-by-Step Execution Flow

1. **Initialization (`Pipeline.\_\_init\_\_`)** from set-up files:

   * Downloads quality-controlled metadata CSVs (`ObservatoryDataFrame`, `OmicsDataFrame`, `SamplingDataFrame`) from the `arms-mbon/data\_workspace` GitHub repository.
   * Downloads the marker- and time-window-specific PEMA OTU/ASV abundance table (`PEMADataFrame`) and, for COI only, the PEMA taxonomic confidence table (`PTAXDataFrame`) from `arms-mbon/analysis\_release\_<N>`.
   * Downloads the corresponding FASTA sequences (`MDAFasta`) from the VLIZ Marine Data Archive (MDA).
   * Initializes three [`Graphoid`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/graphoid.py#L8-L22) output sinks backed by JSON schemas from `data/schemas/`.
2. **Sample-Level Extraction (`Pipeline.\_extract`)**:

   * Identifies all `ReplicateMaterialSampleID` columns in the PEMA abundance table.
   * For each sample ID (`rmsi`), looks up the corresponding row in the **Sampling Event**, **Omics**, and **Observatory** (matched via `ObservatoryID` + `UnitID`) tables.
   * Enforces strict 1-to-1 integrity: if a sample is missing or duplicated in any metadata table, an `AssertionError` is logged and the pipeline moves to the next sample.
   * Filters the PEMA abundance table to retain only OTUs/ASVs with **read count $> 1$** (`self.pema\_df\[self.pema\_df\[self.rmsi] > 1]`).
3. **Occurrence-Level Transformation (`Pipeline.\_etl`)**:

   * Iterates through every present OTU/ASV in the sample and builds a unique `occurrenceID = f"{rmsi}:{OTU}"`.
   * For **COI**, splits the `OTU` identifier (`ASV\_##:amplicon\_hash`) to match the amplicon hash in `PTAXDataFrame` (`tax\_assignments\_\*.tsv`) and filters taxonomic ranks by a minimum confidence threshold (`MIN\_CONFIDENCE = 0.8`).
   * Resolves taxonomy via [`\_walk\_phylogeny`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L61-L105) (detailed below) and queries the ENA Portal API ([`\_get\_ena\_sample\_accession`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L152-L169)) to get from the `sequencing run` accessions the ENA `sample\_accession` IDs (`materialSampleID`).
   * Populates row nodes in `OccurrenceCore`, `DNAExtension`, and `ExtendedMeasurementOrFactsExtension`.
4. **Serialization (`Pipeline.\_serialize`)**:

   * Writes the accumulated Pandas DataFrames inside each `Graphoid` sink to `data/output/<time\_window>\_<genomic\_region>/{occurrence,dnaextension,emof}.csv`.

\---

## 3\. Taxonomic Resolution Logic (`\_walk\_phylogeny`)

Because **EurOBIS** prioritizes marine-verified taxonomic names from **WoRMS** (World Register of Marine Species) or official names in the **Catalogue of Life (CoL)** over raw NCBI taxonomy names, [`\_walk\_phylogeny`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L61-L105) in [function.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py) walks up the taxonomic lineage until it finds an accepted match:

```mermaid
flowchart TD

&#x20;   START(\["Input: PEMA Leaf (Taxon + NCBI ID) \& Lineage Tree<br/>(For COI: lowest rank in PTAX with confidence ≥ 0.8)"])

&#x20;   HAS\_ID{"Has NCBI TaxID?"}

&#x20;   RESOLVE\["Resolve Taxon Name → NCBI TaxID<br/>via ncbi\_taxonomist\_wrapper.resolve()"]

&#x20;   WORMS{"Match in WoRMS API?<br/>(AphiaRecordByExternalID)"}

&#x20;   RET\_WORMS(\["Return WoRMS Aphia as scientificNameID,<br/>scientificName \& taxonRank" + <br/> return 'marine/nonmarine' flags])

&#x20;   COL{"Get info for NCBI ID from ENA Taxonomy +<br/>match name in Catalogue of Life API?"}

&#x20;   RET\_COL(\["Return empty LSID for scientificNameID,<br/>CoL scientificName \& ENA taxonRank"])

&#x20;   HAS\_PARENT{"Parent rank available<br/>from PEMA data?"}

&#x20;   STEP\_UP\["Pop parent taxon from lineage<br/>\& recurse \_walk\_phylogeny() i.e. go up to the preceding rank"]

&#x20;   FALLBACK(\["Fallback:<br/>urn:lsid:marinespecies.org:taxname:1<br/>scientificName = 'Biota'", rank='Kingdom'])



&#x20;   START --> HAS\_ID

&#x20;   HAS\_ID -- Yes --> WORMS

&#x20;   HAS\_ID -- No --> RESOLVE --> WORMS

&#x20;   WORMS -- Yes --> RET\_WORMS

&#x20;   WORMS -- No --> COL

&#x20;   COL -- Yes --> RET\_COL

&#x20;   COL -- No --> HAS\_PARENT

&#x20;   HAS\_PARENT -- Yes --> STEP\_UP --> HAS\_ID

&#x20;   HAS\_PARENT -- No --> FALLBACK```

> \[!NOTE]
> \*\*Special Case Handling in `function.py`\*\*:
> - \*\*Kingdom-only Eukaryota matches\*\*: If the resolved `scientificName` is `"Eukaryota"`, \[`scientific\_name`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms\_toolbox/function.py#L291-L305) reports `"Biota"` with rank `"Kingdom"`.
> - \*\*Terrestrial Taxa Flagging\*\*: \[`occurrence\_remarks`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms\_toolbox/function.py#L254-L268) checks the WoRMS record's `isTerrestrial` flag and appends `"uncertain classification (not marine)"` if set to `1`.



```mermaid

flowchart TD

&#x20;   ENTRY(\["1. Entrypoints (arms\_toolbox/function.py)<br/>• scientific\_name\_id() \[L278]<br/>• scientific\_name() \[L291]<br/>• taxon\_rank() \[L307]"])

&#x20;   COI\_CHECK{"Is ptax present?<br/>(COI marker only)<br/>\[function.py L279/292/308]"}

&#x20;   PTAX\_FILTER\["2. Filter COI Ranks ≥ 0.8 Confidence<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_get\_reliable\_taxon(ptax, min\_confidence=0.8) \[L38]"]

&#x20;   WALK\_START\["3. Start Recursive Phylogeny Walk<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_walk\_phylogeny(pema\_leaf, pema\_tree) \[L61]"]

&#x20;   PARSE\_LEAF\["4. Split 'Taxon:NCBI\_TaxID' String<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_parse\_pema\_leaf(pema\_leaf) \[L16]"]

&#x20;   HAS\_ID{"Has ncbi\_tax\_id?<br/>\[function.py L66]"}

&#x20;   RESOLVE\_NCBI\["5a. Resolve Taxon Name → NCBI TaxID<br/>File: arms\_toolbox/function.py → \_get\_ncbi\_tax\_id(taxon) \[L172]<br/>calls arms\_toolbox/ncbi\_taxonomist\_wrapper.py → resolve(name) \[L15]"]

&#x20;   HAS\_ID\_AFTER{"ncbi\_tax\_id found?<br/>\[function.py L70]"}

&#x20;   WORMS\_CALL\["6. Query WoRMS REST API by NCBI TaxID<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_get\_aphia\_taxonomy(ncbi\_tax\_id) \[L108]"]

&#x20;   WORMS\_MATCH{"Match returned<br/>by WoRMS?<br/>\[function.py L72]"}

&#x20;   RET\_WORMS(\["Return WoRMS Match \[function.py L73-L77]<br/>• scientific\_name\_id = 'urn:lsid:marinespecies.org:taxname:{AphiaID}'<br/>• scientific\_name = atax\['scientificname']<br/>• taxon\_rank = atax\['rank']"])

&#x20;   ENA\_CALL\["7. Query ENA Taxonomy REST API by NCBI TaxID<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_get\_ena\_taxonomy(ncbi\_tax\_id) \[L139]<br/>(Fetches etax\['rank'] \& fallback etax\['scientificName'])"]

&#x20;   COL\_CALL\["8. Query Catalogue of Life (ChecklistBank) API by Name<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_get\_col\_scientific\_name(taxon) \[L121]<br/>(Fetches official col\_scientific\_name)"]

&#x20;   COL\_MATCH{"Does taxon == col\_scientific\_name?<br/>\[function.py L83]"}

&#x20;   RET\_COL(\["Return CoL + ENA Match \[function.py L84-L87]<br/>• scientific\_name\_id = '' (no WoRMS LSID)<br/>• scientific\_name = taxon (validated by CoL)<br/>• taxon\_rank = etax\['rank'] (supplied by ENA)"])

&#x20;   PARSE\_TREE\["9. Extract Ancestor Lineage List<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_parse\_pema\_tree(pema\_tree, taxon) \[L30]<br/>calls \_preprocess\_pema\_tree(classification) \[L24]"]

&#x20;   HAS\_PARENT{"Parent taxon left<br/>in lineage?<br/>\[function.py L93]"}

&#x20;   STEP\_UP\["10. Pop Parent Taxon \& Recurse<br/>File: arms\_toolbox/function.py<br/>Subroutine: \_walk\_phylogeny(taxon=parent, lineage=remaining) \[L98-L100]"]

&#x20;   FALLBACK(\["11. Root Fallback \[function.py L102-L105]<br/>• scientific\_name\_id = 'urn:lsid:marinespecies.org:taxname:1'<br/>• scientific\_name = 'Biota'<br/>• taxon\_rank = ''"])

&#x20;

&#x20;   ENTRY --> COI\_CHECK

&#x20;   COI\_CHECK -- Yes (COI) --> PTAX\_FILTER --> WALK\_START

&#x20;   COI\_CHECK -- No (18S / ITS) --> WALK\_START

&#x20;   WALK\_START --> PARSE\_LEAF --> HAS\_ID

&#x20;   HAS\_ID -- Yes --> WORMS\_CALL

&#x20;   HAS\_ID -- No --> RESOLVE\_NCBI --> HAS\_ID\_AFTER

&#x20;   HAS\_ID\_AFTER -- Yes --> WORMS\_CALL

&#x20;   HAS\_ID\_AFTER -- No --> PARSE\_TREE

&#x20;   WORMS\_CALL --> WORMS\_MATCH

&#x20;   WORMS\_MATCH -- Yes --> RET\_WORMS

&#x20;   WORMS\_MATCH -- No --> ENA\_CALL --> COL\_CALL --> COL\_MATCH

&#x20;   COL\_MATCH -- Yes --> RET\_COL

&#x20;   COL\_MATCH -- No --> PARSE\_TREE

&#x20;   PARSE\_TREE --> HAS\_PARENT

&#x20;   HAS\_PARENT -- Yes --> STEP\_UP --> HAS\_ID

&#x20;   HAS\_PARENT -- No --> FALLBACK

```

> \\\[!Above digram]
> \*\*added the subroutine names to the decision tree\*\*:

&#x20;

\---

## 4\. Complete File-by-File Reference

### Root Directory

|File|Role \& Responsibilities|
|-|-|
|[main.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/main.py)|**Entrypoint script.** Configures logging, sets `DATA\_RELEASE = "002"`, and sequentially runs [`PipelineITS`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_its.py#L3-L5), [`Pipeline18S`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_18s.py#L3-L5), and [`PipelineCOI`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_coi.py#L5-L21) across the configured `time\_window` cohorts (`September2020`, `April2021`, `August2023`) with their direct VLIZ MDA FASTA URLs.|
|[Makefile](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/Makefile)|Convenience task runner providing `make run` (`python main.py`) and `make install` (`pip install -r requirements.txt`).|
|[requirements.txt](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/requirements.txt)|Python dependencies: `pandas`, `biopython`, `pycountry`, `ncbi-taxonomist`, `requests`, `openpyxl`, `tqdm`, and `pyyaml`.|
|[README.md](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/README.md)|Project background, history of Data Releases 001 and 002, taxonomic resolution rationale for OBIS/EurOBIS, and setup instructions for Linux and Windows.|

\---

### Core Package: `arms\_toolbox/`

|File|Key Classes / Functions|Role \& Responsibilities|
|-|-|-|
|[`\_\_init\_\_.py`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/__init__.py)|—|Package marker file for `arms\_toolbox`.|
|[pipeline.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline.py)|[`Pipeline`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline.py#L14-L177)|Base ETL orchestrator. Fetches remote inputs, loops over `ReplicateMaterialSampleID`s and OTUs, invokes field transformation functions, populates the three Darwin Core sinks, and serializes them to CSV.|
|[pipeline\_18s.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_18s.py)|[`Pipeline18S`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_18s.py#L3-L5)|Subclass of `Pipeline` pre-configured with `genomic\_region="18S"`.|
|[pipeline\_its.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_its.py)|[`PipelineITS`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_its.py#L3-L5)|Subclass of `Pipeline` pre-configured with `genomic\_region="ITS"`.|
|[pipeline\_coi.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_coi.py)|[`PipelineCOI`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/pipeline_coi.py#L5-L21)|Subclass of `Pipeline` pre-configured with `genomic\_region="COI"`. Overrides `emof\_properties` with `EMOFPropertiesCOI` and registers three COI-specific eMoF functions (`originalScientificName`, `originalTaxonRank`, `originalScientificNameConfidenceLevel`).|
|[data\_source.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py)|[`ObservatoryDataFrame`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py#L36-L39), [`OmicsDataFrame`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py#L41-L44), [`SamplingDataFrame`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py#L46-L49), [`PEMADataFrame`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py#L51-L72), [`PTAXDataFrame`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py#L74-L99), [`MDAFasta`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_source.py#L102-L115)|Data ingestion adapters. Downloads CSV/Excel/TSV tables from GitHub (`data\_workspace` and `analysis\_release\_<N>`) into Pandas DataFrames and parses FASTA sequences from VLIZ MDA using Biopython (`SeqIO`). Normalizes COI grouped sequence IDs.|
|[data\_sink.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_sink.py)|[`OccurrenceCore`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_sink.py#L4-L8), [`DNAExtension`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_sink.py#L11-L15), [`ExtendedMeasurementOrFactsExtension`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/data_sink.py#L18-L22)|Binds each of the three Darwin Core output tables to its target CSV path (`./data/output/<time\_window>\_<genomic\_region>/\*.csv`) and corresponding JSON schema in `./data/schemas/`.|
|[graphoid.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/graphoid.py)|[`Graphoid`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/graphoid.py#L8-L22)|Schema-backed tabular/graph abstraction. Loads default column values from a JSON schema, appends rows via [`add\_node(\*\*kwargs)`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/graphoid.py#L14-L18), and writes out CSV files via [`serialize()`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/graphoid.py#L20-L22).|
|[function.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py)|[`\_walk\_phylogeny`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L61-L105), [`scientific\_name`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L291-L305), [`dna\_sequence`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L323-L327), [`submerged\_time`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/function.py#L185-L192), etc.|Domain transformation and API integration layer. Houses `@cache`-decorated REST clients for WoRMS, Catalogue of Life, and ENA, plus field-level formatters for countries (ISO alpha-3 to alpha-2), habitats, remarks, PubChem preservative URIs, and DNA sequences.|
|[parameter.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/parameter.py)|`EMOFProperties`, `EMOFPropertiesCOI`|Loads the JSON property templates for standard eMoF measurements and COI-extended eMoF measurements at import time.|
|[ncbi\_taxonomist\_wrapper.py](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/ncbi_taxonomist_wrapper.py)|[`resolve(name)`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/ncbi_taxonomist_wrapper.py#L15-L43)|Programmatic wrapper around the `ncbi-taxonomist` library to resolve a scientific name string to an NCBI TaxID without invoking its CLI.|
|[arms\_toolbox/README.md](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/arms_toolbox/README.md)|—|Technical documentation detailing the `arms\_toolbox` classes, ID formats across markers, and error handling behavior.|

\---

### Schemas \& Outputs: `data/`

|Path|Role \& Responsibilities|
|-|-|
|[data/README.md](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/README.md) \& [data/schemas/README.md](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/README.md)|Overview of the `data/` directory and schema files.|
|[occurrence\_core\_schema\_18s.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/occurrence_core_schema_18s.json), [occurrence\_core\_schema\_coi.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/occurrence_core_schema_coi.json), [occurrence\_core\_schema\_its.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/occurrence_core_schema_its.json)|Define the 34 Darwin Core columns for `occurrence.csv` along with static metadata defaults (e.g., `basisOfRecord: "MaterialSample"`, `institutionCode: "ARMS-MBON"`, marker-specific `identificationRemarks` with classifier/clustering parameters).|
|[dna\_extension\_schema\_18s.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/dna_extension_schema_18s.json), [dna\_extension\_schema\_coi.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/dna_extension_schema_coi.json), [dna\_extension\_schema\_its.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/dna_extension_schema_its.json)|Define the 17 columns for `dnaextension.csv`, including hardcoded primer sequences, primer names, literature references, target gene (`18S`, `COI`, `ITS`), sequencing platform (`MiSeq`), and reference database (`PR2`, `MIDORI`, `UNITE`).|
|[extended\_measurement\_or\_facts\_extension\_schema.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/extended_measurement_or_facts_extension_schema.json)|Defines the 7 columns (`occurrenceID`, `measurementType`, `measurementUnit`, `measurementValue`, `measurementTypeID`, `measurementUnitID`, `measurementValueID`) for `emof.csv`.|
|[extended\_measurement\_or\_facts\_extension\_properties\_schema.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/extended_measurement_or_facts_extension_properties_schema.json)|Lists the 5 standard eMoF measurements (`submergedTime`, `preservativeUsed`, `lowerLimitFilterSize`, `fieldReplicateID`, `technicalReplicateID`) and their BODC NERC Vocabulary URIs (`P01`/`P06`).|
|[extended\_measurement\_or\_facts\_extension\_properties\_schema\_coi.json](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/data/schemas/extended_measurement_or_facts_extension_properties_schema_coi.json)|Extends the 5 standard eMoF measurements with 3 COI-specific measurements: `originalScientificName`, `originalTaxonRank`, and `originalScientificNameConfidenceLevel`.|
|`data/output/<TimeWindow>\_<Marker>/`|Generated output directories (`April2021\_18S`, `April2021\_COI`, `April2021\_ITS`, `August2023\_18S`, `August2023\_COI`, `September2020\_18S`, `September2020\_COI`, `September2020\_ITS`), each containing `occurrence.csv`, `dnaextension.csv`, and `emof.csv`.|

\---

### Post-Hoc \& Historical Scripts: `related\_scripts/`

|File|Role \& Responsibilities|
|-|-|
|[related\_scripts/README.md](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/related_scripts/README.md)|Documents how standalone R scripts and OBIS geographic outlier QC were used to retrospectively fix low-confidence COI assignments and 18S primer references in Data Releases 001 and 002.|
|[ARMS\_issues\_COI\_2018\_2020.R](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/related_scripts/ARMS_issues_COI_2018_2020.R) \& [`.RData`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/related_scripts/ARMS_issues_COI_2018_2020.RData)|Standalone R script and workspace image used to filter low-confidence (< 0.8) COI identifications and relabel Kingdom-only `Eukaryota` hits to `Biota incertae sedis` for **Data Release 001**.|
|[ARMS\_issues\_COI\_2020\_2021.R](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/related_scripts/ARMS_issues_COI_2020_2021.R) \& [`.RData`](file:///C:/Users/cedricd/Documents/Github/darwin-core-generator/related_scripts/ARMS_issues_COI_2020_2021.RData)|Standalone R script and workspace image applying equivalent COI taxonomic confidence corrections for **Data Release 002** (logic now also integrated into `arms\_toolbox/function.py`).|



