
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

Running sniffles with a constant min allele frequency (`--mosaic-af-min 0.05`, the default) shows a peak at 100x instead:

![](figures/26.png)

Calls from the BAM with max coverage:

![](figures/17.png)
![](figures/18.png)




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
* `SMHT001`: cirrotic liver sample sequenced at ~124x with Illumina NovaSeq X. Data downloaded from workspace [SMaHT_Benchmarking_Short_Read](https://app.terra.bio/#workspaces/smaht-gcc-short-read/SMaHT_Benchmarking_Short_Read/data), table `SMAHT001_collaborator_short_read`, field `SMAHT001_collaborator_short_read_id=SM-OLQZF`, fields `collaborator_sample_id=SMHT001-3I-001A1`.

Many short-read somatic callers work in a tumor-normal setting, so we select the following data as "normal":

* `ST001`: a skin sample sequenced at ~100x with Illumina NovaSeq X. Data downloaded from [the benchmarking section of the data portal](https://data.smaht.org/data/benchmarking/donor-st001#skin), file accession [SMAFINU18SZV](https://data.smaht.org/output-files/40ba2ef7-0de0-44d9-992a-189731059bb0/).
* `SMHT001`: a skin sample sequenced at ~90x with Illumina NovaSeq X. Data downloaded from workspace [SMaHT_Benchmarking_Short_Read](https://app.terra.bio/#workspaces/smaht-gcc-short-read/SMaHT_Benchmarking_Short_Read/data), table `SMAHT001_collaborator_short_read`, field `SMAHT001_collaborator_short_read_id=SM-OLRTD`, fields `collaborator_sample_id=SMHT001-3AF-001A2`.

We run [Manta](https://github.com/Illumina/manta?tab=readme-ov-file) in tumor-normal  somatic mode (Manta has also a tumor-only mode, but it is [experimental](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#introduction)). For `ST001` (healthy) we get 134 calls, 120 of which are BNDs. Example:
```
chr1    10012296        MantaBND:1200:0:0:0:0:0:0       G       [chr1:10012368[G        .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:1200:0:0:0:0:0:1;CIPOS=0,3;HOMLEN=3;HOMSEQ=CTT;SOMATIC;SOMATICSCORE=25;BND_DEPTH=83;MATE_BND_DEPTH=85  PR:SR   58,0:88,0       63,2:115,2
chr1    10012365        MantaBND:1200:0:0:0:0:0:1       A       [chr1:10012299[A        .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:1200:0:0:0:0:0:0;CIPOS=0,3;HOMLEN=3;HOMSEQ=GCA;SOMATIC;SOMATICSCORE=25;BND_DEPTH=85;MATE_BND_DEPTH=83  PR:SR   58,0:88,0       63,2:115,2
chr1    16035669        MantaDEL:0:3498:3499:0:0:0      G       <DEL>   .       MinSomaticScore END=16058763;SVTYPE=DEL;SVLEN=-23094;IMPRECISE;CIPOS=-135,136;CIEND=-115,115;SOMATIC;SOMATICSCORE=25    PR      140,1   136,8
chr1    17868650        MantaBND:0:3926:3929:0:0:0:0    A       [chr11:68505958[A       .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:0:3926:3929:0:0:0:1;IMPRECISE;CIPOS=-295,295;SOMATIC;SOMATICSCORE=21;BND_DEPTH=90;MATE_BND_DEPTH=119   PR      105,1   96,9
chr1    21981417        MantaDUP:TANDEM:2443:0:2:0:0:0  C       <DUP:TANDEM>    .       MinSomaticScore END=22006096;SVTYPE=DUP;SVLEN=24679;IMPRECISE;CIPOS=-205,205;CIEND=-152,152;SOMATIC;SOMATICSCORE=10     PR      166,2   135,7
chr1    28377206        MantaBND:3449:0:1:0:0:0:1       A       A]chr1:28377946]        .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:3449:0:1:0:0:0:0;CIPOS=0,18;HOMLEN=18;HOMSEQ=GGATTACAGGCGTGAGCC;SOMATIC;SOMATICSCORE=16;BND_DEPTH=92;MATE_BND_DEPTH=89 PR:SR   101,0:114,0     109,1:139,3
chr1    28377928        MantaBND:3449:0:1:0:0:0:0       T       T]chr1:28377224]        .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:3449:0:1:0:0:0:1;CIPOS=0,18;HOMLEN=18;HOMSEQ=GGCTCACGCCTGTAATCC;SOMATIC;SOMATICSCORE=16;BND_DEPTH=89;MATE_BND_DEPTH=92 PR:SR   101,0:114,0     109,1:139,3
```

For `SMHT001` (cirrotic) we get 32 calls, 26 of which are BNDs. Example:
```
chr1    90378706        MantaBND:8990:0:1:0:0:0:0       T       T]chr7:134929711]       .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:8990:0:1:0:0:0:1;IMPRECISE;CIPOS=-128,128;SOMATIC;SOMATICSCORE=16;BND_DEPTH=12;MATE_BND_DEPTH=53       PR      17,0    28,5
chr2    91518954        MantaDEL:30438:1:1:0:1:0        CTGGCTGCCTGGCTGGCTGTTTGGCTTGGCTGGCTTGGCTGGCTGGGTGGCTTGGCTGGTTTGGCTAGCTGGCTGGCTGGGTAGCTTGGCTAGCTGGCTGGCTTGGCTGGTTGGCTGGCTTGGCTGGCTTGGCTAGCCGGCTGGCTTAGTTGGCTGGATGGCTTGGCTGGCATGCCTGGCT     C       .       MaxMQ0Frac;MinSomaticScore      END=91519134;SVTYPE=DEL;SVLEN=-180;CIGAR=1M180D;CIPOS=0,6;HOMLEN=6;HOMSEQ=TGGCTG;SOMATIC;SOMATICSCORE=14        PR:SR   19,0:58,0       38,0:83,8
chr2    94646936        MantaBND:30805:0:3:0:0:0:0      C       C[chr6:9285873[ .       MinSomaticScore SVTYPE=BND;MATEID=MantaBND:30805:0:3:0:0:0:1;IMPRECISE;CIPOS=-160,161;SOMATIC;SOMATICSCORE=21;BND_DEPTH=72;MATE_BND_DEPTH=41    PR        37,1    32,8
chr4    20885473        MantaBND:62961:0:2:0:0:0:1      G       G]chr5:58943836]        .       PASS    SVTYPE=BND;MATEID=MantaBND:62961:0:2:0:0:0:0;IMPRECISE;CIPOS=-85,85;SOMATIC;SOMATICSCORE=56;BND_DEPTH=103;MATE_BND_DEPTH=113    PR        119,0   107,9
chr4    63202121        MantaDEL:66559:0:0:0:4:0        TATATAGGATATATAGTATATATATAATATATATAGGATATATAGTATATATATAATATAGGATATATAGTATATTTACACTATATATAATGTCTAATAGGATATATAGTATATATATA T       .       MinSomaticScore END=63202239;SVTYPE=DEL;SVLEN=-118;CIGAR=1M118D;CIPOS=0,30;HOMLEN=30;HOMSEQ=ATATAGGATATATAGTATATATATAATATA;SOMATIC;SOMATICSCORE=10        PR:SR   4,0:20,4        7,0:18,5
chr5    58943836        MantaBND:62961:0:2:0:0:0:0      T       T]chr4:20885473]        .       PASS    SVTYPE=BND;MATEID=MantaBND:62961:0:2:0:0:0:1;IMPRECISE;CIPOS=-93,94;SOMATIC;SOMATICSCORE=56;BND_DEPTH=113;MATE_BND_DEPTH=103    PR        119,0   107,9
chr5    77395471        MantaBND:84650:0:4:1:0:0:1      C       C]chr12:95432519]       .       PASS    SVTYPE=BND;MATEID=MantaBND:84650:0:4:1:0:0:0;IMPRECISE;CIPOS=-81,81;SOMATIC;SOMATICSCORE=35;BND_DEPTH=106;MATE_BND_DEPTH=104    PR        100,0   70,7
```

<!--
Other SV callers we should try:
* [GRIDSS](https://github.com/PapenfussLab/gridss?tab=readme-ov-file) in single-sample mode. This is unlikely to outputs somatic calls, since the [readme suggests](https://github.com/PapenfussLab/gridss?tab=readme-ov-file#how-do-i-perform-tumournormal-somatic-variant-calling) to jointly call on all samples from a patient. We should also run `gridss_somatic_filter` at the end, but this seems to [require](https://github.com/PapenfussLab/gridss?tab=readme-ov-file#how-do-i-create-the-panel-of-normals-required-by-gridss_somatic_filter) a panel of normals since it still assumes a tumor/normal setting. 

Delly seems to [require](https://github.com/dellytools/delly#somatic-sv-calling) matched tumor/normal samples, so we don't run it.
-->


## SMHT001 (liver failure)

We IGV'd each one of the top 50 expressed genes in the liver according to [GTEx](https://www.gtexportal.org/home/tissue/Liver), but by eye we could only find intronic SVs in C1R and HPD that are germline (we checked other tissues from the same donor and they appear there as well).

![](figures/27.png)
![](figures/28.png)






---

# Tandem repeat analysis

We IGV'd (in the ST001 230x PacBio BAM) each one of the top 50 expressed genes in the liver according to [GTEx](https://www.gtexportal.org/home/tissue/Liver), and we saw multiple haplotypes in some TRs:

![](figures/29.png)
![](figures/30.png)
![](figures/31.png)
![](figures/32.png)

## Long vs short reads in TRs

It is well known that PacBio reads (top) can capture the structure of TRs better than Illumina reads (bottom). Examples with candidate somatic variation from ST001 230x PacBio:

![](figures/33.png)
![](figures/34.png)
![](figures/35.png)
![](figures/36.png)

![](figures/39.png)
![](figures/37.png)
![](figures/38.png)

![](figures/40.png)
![](figures/41.png)

PacBio reads can also show the 5mC status of TRs (examples from ST001 230x PacBio):

![](figures/42.png)
![](figures/43.png)


## Multiple sequence alignment of the reads

We extract reads from the ST001 230x PacBio BAM with hapestry's command:
```
extract_reads_from_windows --output_dir ./reads_with_q/ \
	--bam_csv ${SAMPLES_LIST} \
	--windows ${FLANKED_WINDOWS_BED} \
	--bam_not_hardclipped \
	--flank_length 0 \
	--fetch_max_length 100000 \
	--get_qualities \
	--tags NM \
	--n_threads 1 \
	--force_forward
```
enabling and disabling `--require_spanning` as needed. Then, we load the FASTA into Jalview and use MAFFT with preset FFT-NS-1 (Speed oriented). Note that a multiple sequence alignment automatically hides all germline homs, and puts all the germline hets in the same cluster.

### chr4:99305502-99305645
We only show reads that span the window:

![](figures/33.png)
![](figures/44.png)

### chr1:16563818-16567204
The following INS seems to be homozygous, since it appears frequently in the BAM and is hidden in the MSA. Only one read seems to have a large INS in the MSA (which shows the right end of the window):

![](figures/39.png)
![](figures/45.png)

### chr19:3971709-3975587
The INS and DEL in the BAM seem to be homozygous, since they do not appear in the MSA. The latter only contains small variation in a homopolymer:

![](figures/35.png)
![](figures/46.png)

### chr19:3988584-3990151
Every alignment in the BAM seems to be affected either by the large INS or by soft-clips. The MSA shows that every spanning read is equivalent. However, non-spanning reads are all suspiciously cut at the same position, which might indicate a fetching error (the white at the top corresponds to reads that align on the right side of the INS but are cut right at the INS).

![](figures/47.png)
![](figures/48.png)


---

## To do

https://pmc.ncbi.nlm.nih.gov/articles/PMC12339610/

https://trtools.readthedocs.io/en/stable/source/prancSTR.html

https://academic.oup.com/bioinformatics/article/40/8/btae485/7723996

TR and liver: https://research.edgehill.ac.uk/en/publications/length-of-variable-numbers-of-tandem-repeats-in-the-carboxyl-este#:~:text=Length%20of%20Variable%20Numbers%20of%20Tandem%20Repeats,M%C3%B6ssner%2C%20Claudia%20Ruffert%2C%20Mario%20Krehan%2C%20Christian%20Zapf

https://pubmed.ncbi.nlm.nih.gov/40004542/

https://www.mdpi.com/2073-4425/16/2/213
