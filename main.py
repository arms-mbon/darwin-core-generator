from arms_toolbox.pipeline_its import PipelineITS
from arms_toolbox.pipeline_18s import Pipeline18S
from arms_toolbox.pipeline_coi import PipelineCOI

# ITS
print("PipelineITSApril2021")
PipelineITS(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c72a7dc04901354082",
).run()

print("PipelineITSJuly2019")
PipelineITS(
    time_window="July2019",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c72a6455c030029062",
).run()

# 18S
print("Pipeline18SApril2021")
Pipeline18S(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c30492304947055810",
).run()

print("Pipeline18SJanuary2020")
Pipeline18S(
    time_window="January2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c304a6079692934739",
).run()

print("Pipeline18SJanuary2022")
Pipeline18S(
    time_window="January2022",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c304974c2824435732",
).run()

print("Pipeline18SJuly2019")
Pipeline18S(
    time_window="July2019",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c3049223a375779320",
).run()

print("Pipeline18SMay2021")
Pipeline18S(
    time_window="May2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c304a04ec570416353",
).run()

print("Pipeline18SSeptember2020")
Pipeline18S(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5c30495306133353334",
).run()

print("Pipeline18SAugust2023")
Pipeline18S(
    time_window="August2023",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_65ddecfdd3c8b943104041",
).run()

print("Pipeline18SARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO")
Pipeline18S(
    time_window="ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_65ddecfdd4483691581806",
).run()

# COI
print("PipelineCOIApril2021")
PipelineCOI(
    time_window="April2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5aec1c61d2304932779",
).run()

print("PipelineCOIJanuary2020")
PipelineCOI(
    time_window="January2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5aec1d37cd080876413",
).run()

print("PipelineCOIJanuary2022")
PipelineCOI(
    time_window="January2022",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5aec1c9cc8783258670",
).run()

print("PipelineCOIJuly2019")
PipelineCOI(
    time_window="July2019",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5aec1c8507131186606",
).run()

print("PipelineCOIMay2021")
PipelineCOI(
    time_window="May2021",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5aec1d22f7293220847",
).run()

print("PipelineCOISeptember2020")
PipelineCOI(
    time_window="September2020",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_64e5aec1ce619790794228",
).run()

print("PipelineCOIAugust2023")
PipelineCOI(
    time_window="August2023",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_65ddecfdd4817905877667",
).run()

print("PipelineCOIARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO")
PipelineCOI(
    time_window="ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO",
    aligned_assignment_url="https://mda.vliz.be/directlink.php?fid=VLIZ_00000615_65ddecfdd3dbd317188875",
).run()

print("process finished")
