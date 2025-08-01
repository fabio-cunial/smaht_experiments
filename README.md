
# Effect of coverage on somatic SV calling

[Metadata of all SMaHT samples](https://docs.google.com/spreadsheets/d/11T_QpVq4XEfupEeGD9IW5oVt1our6u3KzEtP6LBEc5w/edit?usp=sharing) at the time of this study.

# Liver

We consider the following input:
* One healthy liver sample ([donor ST001](https://data.smaht.org/data/benchmarking/donor-st001#liver) in the benchmarking section of the data portal), sequenced at `~92*6=550`x  with PacBio Revio and at `~92*6=550`x with ONT PromethION 24 (coverages estimated from chr1).
* One liver sample from an individual whose cause of death was liver failure and whose death circumstances were alcohol abuse (donor SMHT001 in workspace [SMaHT_Short_Read_Long_Read_Analysis](https://app.terra.bio/#workspaces/smaht-gcc-short-read/SMaHT_Short_Read_Long_Read_Analysis/data), table `SMAHT001_collaborator_long_read`, field `sample_collaborator_id=SMHT001-3I-001A2`), sequenced at ~10x with PacBio Revio.
