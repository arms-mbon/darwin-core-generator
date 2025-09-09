"""
modify the time_window and aligned_assigment_url according to your running of this code
the COI use the all_sequences_grouped fasta files, 18S and ITS use the aligned_assignments ones
"""
from arms_toolbox.pipeline_its import PipelineITS
from arms_toolbox.pipeline_18s import Pipeline18S
from arms_toolbox.pipeline_coi import PipelineCOI

# ITS
print("PipelineITSApril2021")
PipelineITS(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836e1efe1baf206164086"
).run()

print("PipelineITSSeptember2020")
PipelineITS(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836e1efe1b74403284481"
).run()

# 18S
print("Pipeline18SApril2021")
Pipeline18S(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836eb9ad4d41893285473"
).run()

print("Pipeline18SAugust2023")
Pipeline18S(
    time_window="August2023",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836e9574fcf4378071514"
).run()

print("Pipeline18SSeptember2020")
Pipeline18S(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836eb9ae1bae728553767"
).run()

# COI
print("PipelineCOIApril2021")
PipelineCOI(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836ddf0523bd985666107"
).run()

print("PipelineCOIAugust2023")
PipelineCOI(
    time_window="August2023",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836ddf05038e304477491"
).run()

print("PipelineCOISeptember2020")
PipelineCOI(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00001620_6836ec58b6c4d489886842"
).run()
