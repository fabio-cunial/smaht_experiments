
# Effect of coverage on somatic SV calling

[Metadata of all SMaHT samples](https://docs.google.com/spreadsheets/d/11T_QpVq4XEfupEeGD9IW5oVt1our6u3KzEtP6LBEc5w/edit?usp=sharing) at the time of this study.

# Liver

We consider the following input:
* One healthy liver sample ([donor ST001](https://data.smaht.org/data/benchmarking/donor-st001#liver) in the benchmarking section of the data portal), sequenced at ~230x  with PacBio Revio (this includes a ~26x Fiber-seq BAM) and at ~210x with ONT PromethION 24 (coverages estimated from chr1).
* One liver sample from an individual whose cause of death was liver failure and whose death circumstances were alcohol abuse (donor SMHT001 in workspace [SMaHT_Short_Read_Long_Read_Analysis](https://app.terra.bio/#workspaces/smaht-gcc-short-read/SMaHT_Short_Read_Long_Read_Analysis/data), table `SMAHT001_collaborator_long_read`, field `sample_collaborator_id=SMHT001-3I-001A2`), sequenced at ~10x with PacBio Revio.

For each ST001 technology, we merge all the BAMs, we take random samples at multiples of 10x, and on every such subsample we run `sniffles --mosaic` requiring just two reads to support a call (in particular, we set `--mosaic-af-min` as a function of coverage). We only output calls that Sniffles considers somatic, i.e. that have `--mosaic-af-max=0.218` (the default). Each call is annotated with the IDs of the reads that support it.

## SMHT001
There are only 14 calls, all of which except one occur in a TR. 

### Potential candidates

![](figures/5.png)
![](figures/6.png)
![](figures/8.png)
![](figures/10.png)
![](figures/14.png)

### Calls unlikely to be somatic

![](figures/1.png)
![](figures/2.png)
![](figures/9.png)

### Calls near complex events

![](figures/3.png)
![](figures/7.png)
