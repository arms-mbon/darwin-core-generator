"""
# ITS and 18S
Read in the ARMS-MBON data input -- PEMA files (OTU tables), setup files (for the emof), and observatory, omics, and sampling event 
 files -- and turn those into DwCA (occurrence, emof, and dnaextension CSV files).
The emof schema can be optimised to the gene type.
This code is designed to be run over several input files whoes name differ by gene type (18S and ITS) and date (e.g. July2020).
The output is added to data/outputs.
Checks for duplicates and missing information are performed - for duplicates, a warning is logged and the first only is taken.
A check on ENA to find the sample accession numbers for the run accession numbers that we have in the omics input is done.
TODO
I get reports for 18S that there are run accession numbers with no sample accession numbers found in ENA
in fact that is not the case, it misses only some of the occurrence IDs
We only look for WoRMS AphiaIds for the scientitic name level in the data here - we do not travel along the taxon tree to find a match at a higher and higher level. 
To be considered if this should be done - would result in more matches to WoRMS

# COI
Read in the ARMS-MBON data input -- PEMA files (OTU tables), setup files (for the emof), and observatory, omics, and sampling event 
 files -- and turn those into DwCA (occurrence, emof, and dnaextension CSV files).

To be improved
- Having the NCBI Id - APhia ID search put results in a dictionary so that they don't need to be done each time (? unsure if already done?)
- Reading from the tax_assigments files for COI: all our results go to species level, but this may not be so need to take the scientific name, 
  confidence level, and rank, from the entry that is the most specific for that row (not assume that it will be always the last row)
- There are now so many IF COI ELSE in here, that it needs to be split into two scripts (as was originally the case)
- I get reports for 18S that there are run accession numbers with no sample accession numbers found in ENA
  in fact that is not the case, it misses only some of the occurrence IDs
- We only look for WoRMS AphiaIds for the scientitic name level in the data here - we do not travel along the taxon tree to find a match at a higher and higher level. 
  To be considered if this should be done - would result in more matches to WoRMS

#
modify the time_window and aligned_assigment_url according to your running of this code
the COI use the all_sequences_grouped fasta files, 18S and ITS use the aligned_assignments ones
"""
import logging
from arms_toolbox.pipeline_its import PipelineITS
from arms_toolbox.pipeline_18s import Pipeline18S
from arms_toolbox.pipeline_coi import PipelineCOI

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

DATA_RELEASE = "002"
DATASET_NAME = f"data_release_{DATA_RELEASE}"

# ITS
PipelineITS(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836e1efe1baf206164086",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

PipelineITS(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836e1efe1b74403284481",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

# 18S
Pipeline18S(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836eb9ad4d41893285473",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

Pipeline18S(
    time_window="August2023",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836e9574fcf4378071514",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

Pipeline18S(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836eb9ae1bae728553767",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

# COI
PipelineCOI(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836ddf0523bd985666107",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

PipelineCOI(
    time_window="August2023",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836ddf05038e304477491",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()

PipelineCOI(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836ec58b6c4d489886842",
    dataset_name=DATASET_NAME,
    data_release=DATA_RELEASE,
).run()
