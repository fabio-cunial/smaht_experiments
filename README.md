
# Effect of coverage on somatic SV calling

[Metadata of all SMaHT samples](https://docs.google.com/spreadsheets/d/11T_QpVq4XEfupEeGD9IW5oVt1our6u3KzEtP6LBEc5w/edit?usp=sharing) at the time of this study.

# Liver, long reads.

We consider the following input:

* `ST001`: healthy liver sample, sequenced at ~230x with PacBio Revio (this includes a ~26x Fiber-seq BAM) and at ~210x with ONT PromethION 24 (coverages estimated from chr1). Data downloaded from [the benchmarking section of the data portal](https://data.smaht.org/data/benchmarking/donor-st001#liver).
* `SMHT001`: death caused by liver failure, alcohol abuse as a death circumstance. Sequenced at ~10x with PacBio Revio. Data downloaded from workspace [SMaHT_Short_Read_Long_Read_Analysis](https://app.terra.bio/#workspaces/smaht-gcc-short-read/SMaHT_Short_Read_Long_Read_Analysis/data), table `SMAHT001_collaborator_long_read`, field `sample_collaborator_id=SMHT001-3I-001A2`.

For each ST001 technology, we merge all the BAMs, we take random samples at multiples of 10x, and on every such subsample we run [`sniffles --mosaic`](https://github.com/fritzsedlazeck/Sniffles#c-mosaic-sv-calling-non-germline-or-somatic-svs) **requiring just two reads** to support a call (i.e. we set `--mosaic-af-min` as a function of coverage). We only output calls that Sniffles considers as somatic, i.e. that have `--mosaic-af-max=0.218` (the default). Each call is annotated with the IDs of the reads that support it. We only consider calls in the standard chromosomes.

![](figures/15.png)
![](figures/16.png)
Circles: 10x SMHT001 (liver failure). In the 230x ST001 VCF, TR calls do not seem to be enriched in specific TR intervals (according to `bedtools intersect -c`).

<img src="./figures/25.png" width="400" height="400"/>

![](figures/17.png)
![](figures/18.png)

Running sniffles with a constant min allele frequency (`--mosaic-af-min 0.05`, the default) shows a peak at 100x instead:

![](figures/26.png)


## ST001 PacBio (healthy)
Examples of calls of length >6k in the 230x VCF:

![](figures/19.png)

![](figures/20.png)
![](figures/21.png)
![](figures/22.png)
![](figures/23.png)
![](figures/24.png)



## SMHT001 PacBio (liver failure)
There are only 14 total calls, all of which except one occur in a TR. 

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




# Liver, short reads.

We consider the following input:

* `ST001`: healthy liver sample, sequenced at ~122x with Illumina NovaSeq X. Data downloaded from [the benchmarking section of the data portal](https://data.smaht.org/data/benchmarking/donor-st001#liver), file accession [SMAFILE7Y4Y9](https://data.smaht.org/output-files/b784cb52-f497-4113-b341-813ee0e6d700/).
* `SMHT001`: Sequenced at ~124x with Illumina NovaSeq X. Data downloaded from workspace [SMaHT_Benchmarking_Short_Read](https://app.terra.bio/#workspaces/smaht-gcc-short-read/SMaHT_Benchmarking_Short_Read/data), table `SMAHT001_collaborator_short_read`, field `SMAHT001_collaborator_short_read_id=SM-OLQZF`, fields `collaborator_sample_id=SMHT001-3I-001A1`.

We run the following SV callers:
* [Manta](https://github.com/Illumina/manta?tab=readme-ov-file) in single-BAM somatic mode ([experimental](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#introduction)).
* [GRIDSS](https://github.com/PapenfussLab/gridss?tab=readme-ov-file) in single-sample mode. This is unlikely to outputs somatic calls, since the [readme suggests](https://github.com/PapenfussLab/gridss?tab=readme-ov-file#how-do-i-perform-tumournormal-somatic-variant-calling) to jointly call on all samples from a patient. We should also run `gridss_somatic_filter` at the end, but this seems to [require](https://github.com/PapenfussLab/gridss?tab=readme-ov-file#how-do-i-create-the-panel-of-normals-required-by-gridss_somatic_filter) a panel of normals since it still assumes a tumor/normal setting. 

Delly seems to [require](https://github.com/dellytools/delly#somatic-sv-calling) matched tumor/normal samples, so we don't run it.



## SMHT001 (liver failure)

I IGV'd each one of the top 50 expressed genes in the liver according to [GTEx](https://www.gtexportal.org/home/tissue/Liver), but I could only find SVs in C1R and HPD that are germline (I checked other tissues from the same donor and they appear there as well).

![](figures/27.png)
![](figures/28.png)
