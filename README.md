
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




# Note on clipped alignments

In the ST001 230x PacBio BAM, ~3% of all reads have at least one clipped alignment, defined as an alignment with either a 100bp left or a 100bp right soft clip and that involves a canonical chromosome. Often these reads map to more than one chromosome, or contain sequence that does not map to any region of GRCh38. We assume that such reads are chimeras and discard them, unless their clips are aligned and their patterns of mismatching bases are supported by more than one read. The same phenomenon occurs in PacBio samples from HPRC, and in ONT reads from ST001. 

More info on chimeras: Heinz et al. [Detecting foldback artifacts in long reads](https://www.biorxiv.org/content/10.1101/2025.07.15.664946v2.abstract). bioRxiv 2025.




# Tandem repeat analysis

People in the SMaHT network are already focusing on TRs, see e.g. this [concept sheet](https://docs.google.com/document/d/1uVWeXBupx8FzR1l9lGAjPxttVp2uE_O-UgKcprPHlZM/edit?tab=t.0).

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




## Detailed analysis of some TR regions

In the ST001 230x PacBio BAM. 

**Note: some interesting regions in or around these genes might be missing from what follows. Go over IGV in detail once again.**


### ⭐️ ADH1B: TR in the 5' UTR.
<!-- chr4:99305502-99305645 -->

![](figures/33.png)

We extract spanning reads from the ST001 230x PacBio BAM with hapestry's command:
```
extract_reads_from_windows --output_dir ./reads_with_q/ \
	--bam_csv ${SAMPLES_LIST} \
	--windows ${FLANKED_WINDOWS_BED} \
	--bam_not_hardclipped \
	--require_spanning \
	--flank_length 0 \
	--fetch_max_length 100000 \
	--get_qualities \
	--tags NM \
	--n_threads 1 \
	--force_forward
```
Then, we load the FASTA into Jalview and use MAFFT with preset FFT-NS-1 (Speed oriented). Note that a multiple sequence alignment automatically hides all germline homs, and puts all the germline hets in the same cluster.

![](figures/44.png)

This TR has no 5mC marks.

TRGT detects two haplotypes in one interval:
```
chr4	99305501	.	

CATATATATATATATATATATATACAATCACTTAACTATATGAACCACTTGCCCCATAGTGTATATATATATATATATATATATATATA	

CATATATATATATATATATATATATACAATCACTTAACTATATGAACCACTTGCCCCATAGTGTATATATATATATATATATA,
CATATATATATATATATATATATACAATCACTTAACTATATGAACCACTTGCCCCATAGTGTATATATATATATATATATATATATATATA	

.	.	

TRID=4-99305501-99305524-AT,4-99305558-99305562-GT;END=99305589;MOTIFS=AT,GT;STRUC=<VC181780>	
GT:AL:ALLR:SD:MC:MS:AP:AM	
1/2:82,90:66-86,56-95:114,136:25_2,29_2:0(0-24)_0(27-29)_0(38-42)_0(56-58)_1(59-63)_0(63-81),0(0-22)_0(25-27)_0(36-40)_0(54-56)_1(57-61)_0(61-89):0.682927,0.711111:.,.
```
And 250 reads span them (coverage in the region: 330). However, aligning the reads to the two haplotypes shows recurrent CIGAR patterns:

![](figures/115.png)
![](figures/116.png)

There is an overlapping TRGT catalog interval in this region, which is also genotyped by TRGT with two haplotypes:

```
chr4	99305561	.	

GTATATATATATATATATATATATATATATATATATATATATATATATA	

GTATATATATATATATATATATATATATATATATATATATA,
GTATATATATATATATATATATATATATATATATATATATATATATATATA	

.	.	
TRID=4-99305561-99305609-TA;END=99305609;MOTIFS=TA;STRUC=<TR2721601>	
GT:AL:ALLR:SD:MC:MS:AP:AM	
1/2:40,50:24-44,16-54:113,137:20,25:0(0-40),0(0-50):1.000000,1.000000:.,.
```
And 250 reads span them. However, as above, aligning the reads to the two haplotypes shows recurrent CIGAR patterns.

Sniffles: no support.


### ❓NBPF1: TR overlapping exons.
<!-- chr1:16563818-16567204 -->

![](figures/49.png)

This TR seems to be methylated, and to have >=2 distinct 5mC patterns (gray: FWD strand, green: REV strand):

![](figures/52.png)

Zooming in, several clipped alignments seem to have mismatching bases aligned with a hom INS, and the patterns of mismatch seem to be supported by more than one read. These may just be the beginning/end of the same hom INS, but it is worth investigating.

![](figures/39.png)

We extract all alignments (spanning and non-spanning) with the "Export alignments" feature of IGV. Then, we build a POA graph using abPOA in local alignment mode with flags `-m 1 --amb-strand --sort-by-len --result 3` (the TR sequence in the reference is in red):

![](figures/110.png)
<!-- ![](figures/55.png) -->
<!-- ![](figures/56.png) -->

Strangely the reference sequence is not aligned with the INS, the left INS is supported by ~10 reads and the right INS is supported by ~5 reads. Probably the POA graph is inaccurate in this region.

TRGT detects two haplotypes:
```
chr1	16565544	.	

CAGCTCAGTAATGGCCACTTGGAGCAGGAATATGATCTTTATATGGAAGACTCAGTGGATCCTTATCACCTTCATAGAAAGGTACTCACCTCCCACGTCAAGAGAAAAGCCAACATGTTTTTCCTCCAATGCATAAAAGGAACTTCCATAGGGCAGGCAGGAGTCAGGCTGTTCAAGACAACTGGAAGGAGTTGAATAACATCTATCCAGTGAGTCCTGCAAGACTTCAGGCTCTACTGCCTCCAGCAGCTCCCTGCTGAGCCTGGAAAAGTAGGAAAAAGTAAAGAATAAGCCAGGGGGAATCAGAAACCACACAGCCCCAGCTACATTTCATGGCTAACATAAGGAACTGTTTAAACAGAAAAAGGACAGATCCATTAATGAGGTAATGAATTATTGCCTTTATGTTGGGATAGACCAGGGCCAGGTAGAAAAGAATGAAAGAGAAAGACAGGGAGAGGGAGAGGGAGAGAGAGACAGAGGAGAAAGTGAGCTCAGCGAATTGGCCGGGTGACACACTGACGAAGGGGTCAAAGGACACTCTGAGTTAGTGCCCTCGGGACACACAGAGAACAGTGATCATGAAAAGAGTGGGCTCAATAATTTTCCATAAACTTGCTTAAGATTCCATGCAGTTGCCATACAGCCTTTGAGGTATGGTCAACCTACAGTAAGTTAGTAAATGATAAGGGGAGGAAGAAATGGAAACCTAAACATCTACTGCAAGGAAAACCAACAGCAATGTCAGTAGGAGTAATTCAACCTTCGTTGAAAACATGAAATTGAACATACTCTTGTTTTCCCTGGACCTGGCATCTCCAGGTGTCAACACAGAATTAAGCATCCATAATTGCTCAAAGTTACCTGGGGCATGATGGGTCTTGGTCTTCTTCCACTTCTTGGTACTTTTCAATTTCTGCAATAAGTTCAGACATGGACAGACATATTAAGCTGGTTCTCCTACACACATAACAATCCACTGTCTAATCCTCACGCAGGGACTTCAGGCTCCTCAGCATGAGAATAGGACACTGTGAGAGATCTTCTTCAGGAGGCCTGAAGGCTGATCATGATAGAGATTCCTGGGTTTTTGTCCCAGAAACTGTGGGTAAAATTCCCTATTCTGGTAGATCGTTATCCCAAGATCATTTGTCCCAAGTTTGTGCAAATGGTTATGCCATATTTTTCCAATCGATTTAAAGCAAATGCCCCCAAATGGTTGCTGGGAGAAAAACTGCAATATTCAGCCCTGTCTCATCAAATACTCAGATTCTTCATGGTAGCGAGGATTTTAGATGCTGAAATTAGAGTGAAGGATGAAATCTACAAGATCTACAAAATTGAGACAAAATCAGAGTTGTGTGAATTTGTCACATCTGCCCAGATCCAACATCTTGAGAGTGGGATTAGGGTGCCACAGGCATGGCCTGAGACTAGGAAGAGAGCCCTGCTCACTGACCCATCCCTTGCCTGGGCTTCCAAGTGGAACTAGAGTTTCATTCAACCTACATGTGCCTATAGGTCCTCCCTGTGGCAATGACATCTCTCAGCTCAGTAAGGGCCATTTGCAGTAGGAATATGACCCTAACCAGAAGACTCAGTGGATCCTTATCACCTTCATAGAAAGGTACTCACCATCCATGTCAAGAGCCCAGCCAACACGCTGTTGCTCCAATATGTAAAAGGCACTTCTGTAGGGCTGGCATGAGTCAGTCAGTTCAAGATAACCTGAAGGAGTTGAATAACATCTATCCAGTGAGTCCTGCAAGACTTCAGGCCCTTTCTCATCCAGCAGCTCCCTGCTGAGCCTGGAACAGTGGGAAAAAGTAAAGAATAAGCCAGGGGGAATCAGAAACCACACAGCCCCAGCTAGATTTCATGGCTAACATAAGGAAGAGTTTGAAAAGAAAAAGGACAGATCCATTAATGAGGTAACAAATTATTGCCTTTATATTGGGATAGACTAGGGCCAGGTAGAAAAGGATGAAAGAGAAAGACACACACACACACACACACACACACACACACACACACACACACACAGAGTGAGCTCAGTGAATTGGCCAGGTGACACACTGATGAGGGAGTCAACGGTCATTCTCTATTTGTGCTCTCAGGACACACAGTGAACAGTGATCATGAAAAGCATGGCCTCAATAATTTTGCATAAAATGTGCTCAAGTTTCCCTGCAGCCACCATGAGAATACAGCTTTTGAGGTATGGTCAACCTTCACTAGGTTAGTAAATGATAAGGGTAGGAAGAAATGGAAACCTAAACATTTACTCTAATGAGAACCAAAA	

CAGCTCAGTAAGGGCCATTTGCAGTAGGAATATGACCCTAACCAGAAGACTCAGTGGATCCTTATCACCTTCATAGAAAGGTACTCACCATCCATGTCAAGAGCCCAGCCAACACGCTGTTGCTCCAATATGTAAAAGGCACTTCTGTAGGGCTGGCATGAGTCAGTCAGTTCAAGATAACCTGAAGGAGTTGAATAACATCTATCCAGTGAGTCCTGCAAGACTTCAGGCCCTTTCTCATCCAGCAGCTCCCTGCTGAGCCTGGAACAGTGGGAAAAAGTAAAGAATAAGCCAGGGGGAATCAGAAACCACACAGCCCCAGCTAGATTTCATGGCTAACATAAGGAAGAGTTTGAAAAGAAAAAGGACAGATCCATTAATGAGGTAACAAATTATTGCCTTTATGTTGGGATAGACTAGGGCCAGGTAGAAAAGGATGAAAGAGAAAGACACACACACACACACACACACACACACACACACACACACACACAGAGTGAGCTCAGTGAATTGGCCAGGTGACACACTGATGAGGGAGTCAACGGTCATTCTCTATTTGTGCTCTCAGGACACACAGTGAACAGTGATCATGAAAAGCATGGCCTCAATAATTTTGCATAAAATGTGCTCAAGTTTCCCTGCAGCCACCATGAGAATACAGCTTTTGAGGTATGGTCAACCTTCACTAGGTTAGTAAATGATAAGGGTAGGAAGAAATGGAAACCTAAACATTTACTCTAATGAGAACCAAAA,

CAGCTCAGTAATGGCCACTTGGAGCAGGAATATGATCTTTATATGGAAGACTCAGTGGATCCTTATCACCTTCATAGAAAGGTACTCACCTCCCACGTCAAGAGAAAAGCCAACATGTTTTTCCTCCAATGCATAAAAGGAACTTCCATAGGGCAGGCAGGAGTCAGGCTGTTCAAGACAACTGGAAGGAGTTGAATAACATCTATCCAGTGAGTCCTGCAAGACTTCAGGCTCTACTGCCTCCAGCAGCTCCCTGCTGAGCCTGGAAAAGGAGGAAAAAGTAAAGAATAAGCCAGGGGAAATCAGACACAACAGAGCCCCAACTAGGTTTCATGGGTAGCATAAGGAAGTGGTTAAAAAAGTAAAAGGATAGATCCATTAATGAGGTAACAAATTATTGCCTTCATGTTGGGACAGAACAGGGCCAAATGGAAAAGAATGAAAGAGAAAGACAGATAGACACACACACACACACACACACACACACACAGAGAGAGAGAGAATGAGCTCAGTGAATTGTCCAGGTGACACACTGATGAGGGAGTAACAGGACACTCTGAGTTAGTGCCCTCAGGACACACAGCATACAGTGATCAGGAAAGGACTGTGCTCAATAATTTTCCATAAAATGTGCTCAAGTTTCCATGCAGTCGCCATGAGAATACAGTTTTTGAAGTCTGGTCCACCTACAGTAGGTTAGTAAATGATAAGGGGAGGAAGAAATGGAAACCTAAATATCTACTGCAATGAAAACCAACAGCAATGTTAGTAGGAATAATTCAGGCTCGGTTGAAAAGATGTAATCGATAATGTCAGCCCGCCCTGTTTTCCCTGAACCAGGAGTCTCCAGATGTCAACACAGAAGAAGCTGTTCACAATTGCTCAGTTACCTGGGGCATGGTGGGCCTTGGTCTTCTTCCTCTTCTTGGTCCTTTTTAATTCCTGCAATACATTCAGACAGGGACAGACAAAATAAGCCAATTCACCTACACCCGTAACAGTCCACTGTCTAATCCCCACACAGGGATCTCAGGCTCCTCAGCAAGAGAACAGGACAATGTGAGAGATATACTTCAGGAGGCCTGAAAGCTGGTCATGATATTCTTTGGTTTGCATCTCAGAACCAAGGGTGAAATATCCCTATTCTGGTAGATCGTTATCCCAAAATCATTTATCCCAAGTTTGTGCAAACAGTTATGCCTTATTGTTCCCATCAGTTCAAAGACAATGCCCCAGATGATTTCTAGGAGGAAAACTGCAGTATTCAGCCCTGTCTCATCAAATGCCCAGCTCGTTCATGGATGCAAGAATTTTAGACACTGAAATTAGAATGAGGGAGGAAATCTACAAACCCTTGAGTCCAAATCATAGTTCTGTGAATTTTTTACATCTGCCTGGGTCCAATGTGCTGAGAGCGGGCTCAGGTTGATACAGGCATGGCTGGAGACTAGGAATAGAGCCTTGCTCACTGACCCATTTCATGTCTAGGCTTCCAGCGGAGACTACAGTTTCATTACAACCTATATGCGCCCATAGGTCCTGCCTGCGGCAATGACATCTCTCGGGTCAGTAAGGGCCACTGGGAACAGGAATATCACCCCTATCTGGAAGACCAGGTGGAGGCTTATCACCTTCATAGTAAGGTACTCACTGTCCACGTCAAGAGCCAAGCCAAGGTACTGTTCCTCCAATGAGTAAACAGCACTTCTGTAGGGCTGGCCTAAGTCAGGCAGTTCAAGATAACCTGAAGGAGTCGAATAACATCTATCCAGTGAGTCCTGCAAGACTTCAGGCTCTTTCTCATCCAGCAGCTCCCTGCTGAGCCTGGAAAAGTAGGAAAAAGTAAAGAATAAGCCAGGGGGAATCAGAAACCACACAGCCCCAGCTACATTTCATGGCTAACATAAGGAACTGTTTAAACAGAAAAAGGACAGATCCATTAATGAGGTAATGAATTATTGCCTTTATGTTGGGATAGACCAGGGCCAGGTAGAAAAGAATGAAAGAGAAAGACAGGGAGAGGGAGAGGGAGAGAGAGACAGAGGAGAAAGTGAGCTCAGCGAATTGGCCGGGTGACACACTGATGAAAGGGTCAAAGGACACTCTGAGTTAGTGCCCTCGGGACACACAGAGAACAGTGATCATGAAAAGAGTGGGCTCAATAATTTTCCATAAACTTGCTTAAGATTCCATGCAGTTGCCATACAGCCTTTGAGGTATGGTCAACCTACAGTAAGTTAGTAAATGATAAGGGGAGGAAGAAATGGAAACCTAAACATCTACTGCAAGGAAAACCAACAGCAATGTCAGTAGGAGTAATTCAACCTTCGTTGAAAACATGAAATTGAACATACTCTTGTTTTCCCTGGACCTGGCATCTCCAGGTGTCAACACAGAATTAAGCATCCATAATTGCTCAAAGTTACCTGGGGCATGATGGGTCTTGGTCTTCTTCCACTTCTTGGTACTTTTCAATTTCTGCAATAAGTTCAGACATGGACAGACATATTAAGCTGGTTCTCCTACACACATAACAATCCACTGTCTAATCCTCACACAGGGACTTCAGGCTCCTCAGCATGAGAATAGGACACTGTGAGAGATCTTCTTCAGGAGGCCTGAAGGCTGATCATGATAGAGATTCCTGGGTTTTTGTCCCAGAAACTGTGGGTAAAATTCCCTATTCTGGTAGATCGTTATCCCAAGATCATTTGTCCCAAGTTTGTGCAAATGGTTATGCCATATTTTTCCAATCGATTTAAAGCAAATGCCCCCAAATGGTTGCTGGGAGAAAAACTGCAATATTCAGCCCTGTCTCATCAAATACTCAGATTCTTCATGGTAGCGAGGATTTTAGACGCTGAAATTAGAGTGAAGGATGAAATCTACAAGATCTACAAAATTGAGACAAAATCAGAGTTGTGTGAATTTGTCACATCTGCCCAGATCCAACATCTTGAGAGTAGGATTAGGGTGCCACAGGCATGGCCTGAGACTAGGAAGAGAGCCCTGCTCACTGACCCATCCCTTGCCTGGGCTTCCAAGTGGAACTAGAGTTTCATTCAACCTACATGTGCCTATAGGTCCTCCCTGTGGCAATGACATCTCTCAGCTCAGTAAGGGCCATTTGCAGTAGGAATATGACCCTAACCAGAAGACTCAGTGGATCCTTATCACCTTCATAGAAAGGTACTCACCATCCATGTCAAGAGCCCAGCCAACACGCTGTTGCTCCAATATGTAAAAGGCACTTCTGTAGGGCTGGCATGAGTCAGTCAGTTCAAGATAACCTGAAGGAGTTGAATAACATCTATCCAGTGAGTCCTGCAAGACTTCAGGCCCTTTCTCATCCAGCAGCTCCCTGCTGAGCCTGGAACAGTGGGAAAAAGTAAAGAATAAGCCAGGGGGAATCAGAAACCACACAGCCCCAGCTAGATTTCATGGCTAACATAAGGAAGAGTTTGAAAAGAAAAAGGACAGATCCATTAATGAGGTAACAAATTATTGCCTTTATGTTGGGATAGACTAGGGCCAGGTAGAAAAGGATGAAAGAGAAAGACACACACACACACACACACACACACACACACACACACACACACAGAGTGAGCTCAGTGAATTGGCCAGGTGACACACTGATGAGGGAGTCAACGGTCATTCTCTATTTGTGCTCTCAGGACACACAGTGAACAGTGATCATGAAAAGCATGGCCTCAATAATTTTGCATAAAATGTGCTCAAGTTTCCCTGCAGCCACCATGAGAATACAGCTTTTGAGGTATGGTCAACCTTCACTAGGTTAGTAAATGACAAGGGTAGGAAGAAATGGAAACCTAAACATTTACTCTAATGAGAACCAAAA
```
And 219 reads span them (coverage in the region: 340). However, aligning the reads to the two haplotypes shows recurrent CIGAR patterns:

![](figures/117.png)
![](figures/118.png)

Sniffles: no support.


### ❓EEF2
<!-- chr19:3971709-3975587 -->

The ends of clipped alignments seem to be consistently located, and the mismatching bases seem to be supported by more than one read:

![](figures/50.png)

Moreover, this TR seems to be methylated with a consistent pattern, which seems to diverge in the clipped alignments (gray: FWD strand, green: REV strand):

![](figures/53.png)

We extract all alignments (spanning and non-spanning) with the "Export alignments" feature of IGV. Then, we build a POA graph using abPOA with flags `--amb-strand --sort-by-len --result 3` (the TR sequence in the reference is in red):

![](figures/108.png)
<!-- ![](figures/57.png) -->
<!-- ![](figures/58.png) -->

No read supports the red branching arm (only the ref supports it), and the blue branching arm is supported by ~150 reads: this seems to contradict the BAM.

TRGT detects only the hom DEL:
```
chr19	3973097	.	

TGCCTGTAGTCCCGGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCAAGATTGAGCCACTGCACTCCAGCCTGGGTGACAGAGCAAAGACTGTCTCAAAAAAAAAAAAAAATTCAGCTTCCACATACCATCTTAACGTTTGCACTACTAAAAAAGAATAGGACCCGGTATGGTGGCTCACACCTGTAACCCCAGCACTTTGGGAGGCTGAGGCGGGCTGATCTCTTACGGTCGAGCTTGAGACCAGCCTGGCCAACATGGTGAAACCCCATCTCACCTAAATATACAAAAATTACACAGACCTGGTGATAGGCGCCTGTAATTCCAGCTACTCGGGAGGCTGAAGCACAAGAATTGCTTGAGCACTTGAATCCGGAGGCAGAGGCTATAGTGAGCCAAGATCACAACACTGCACTCCATCCTGGGCGATAGAGTGAGATTCTCGAAACACAGTAATCCGTTACCCTGACTGCTGGAAACAGCGTCAGTTACAGCCCCAATTCAAATGCGACAGCCAAACCAAACCTCAGGCATACCTACATGTGGGCAAAACCTCAGCTGAGGCTTGGCCTTGACAATGCAGGGCCAGAGCTGCCCCCGCCTTTGCCGAGAACCCCACAGCCCTACTTCACAGTCCCCACCTGCCTTCCTGAGGACAGAGTAAAAGGCGGAAAATGATCATGTTAACCTTATCAGCACTGACAAGATTAGGCCTCCTCACATCTTAGGCACTATTTAGAAAAACAAATAGTTTTCTAAGCTATTGGCGCTCACACCTGAAATCCCAGCACTTTGGGAGGCTGAGGCAGAAGGATCGCTTGAGCCCGGGAAGGAGACCAGCCTGAGCAACATAGTGAGACCCCATTACTACAAAAAAAGTAAAAATTAGGTCAGCACAGTGGCTCACGCCTATAATTCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACGAGGTCAGGAGATCGAGACTATCCTGGCTAACACGGTGAAACCCTGTCTCTACTAATAAATACAAAAAATTAGCCGGGCATGGTGGCGGGCGCCTGTAGTCCCGGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCAAGATTGAGCCACTGCACTCCAGCCTGGGTGACAGAGCAAAGACTGTCTCAAAAAAAAAAAAAAAT	

TGCCTGTAGTCCCGGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCAAGATTGAGCCACTGCACTCCAGCCTGGGTGACAGAGCAAAGACTGTCTCAAAAAAAAAAAAAAAT	

.	.	

TRID=19-3973223-3973238-A,19-3974295-3974310-A;END=3974311;MOTIFS=A;STRUC=<VC111219>	
GT:AL:ALLR:SD:MC:MS:AP:AM	
1/1:142,142:141-144,141-143:80,66:43,43:0(6-7)_0(16-17)_0(23-24)_0(29-30)_0(33-34)_0(36-37)_0(38-40)_0(47-49)_0(55-56)_0(61-62)_0(68-69)_0(72-73)_0(76-78)_0(79-80)_0(83-84)_0(87-88)_0(92-93)_0(97-98)_0(107-108)_0(109-110)_0(111-112)_0(114-117)_0(118-119)_0(126-141),0(6-7)_0(16-17)_0(23-24)_0(29-30)_0(33-34)_0(36-37)_0(38-40)_0(47-49)_0(55-56)_0(61-62)_0(68-69)_0(72-73)_0(76-78)_0(79-80)_0(83-84)_0(87-88)_0(92-93)_0(97-98)_0(107-108)_0(109-110)_0(111-112)_0(114-117)_0(118-119)_0(126-141):0.302817,0.302817:0.59,0.59
```
And 146 reads span it (coverage in the region: 160). However, aligning the reads to the two haplotypes shows recurrent CIGAR patterns, e.g.:

![](figures/119.png)

Sniffles: no support.


### chr19:3988584-3990151

Several clipped alignments seem to have mismatching bases aligned with a het INS, and the patterns of mismatch seem to be supported by more than one read. These may just be the beginning/end of the same het INS, but it is worth investigating.

![](figures/51.png)

Moreover, this TR seems to be methylated with a consistent pattern, which seems to diverge in the clipped alignments (gray: FWD strand, green: REV strand):

![](figures/54.png)

We extract all alignments (spanning and non-spanning) with the "Export alignments" feature of IGV. Then, we build a POA graph using abPOA with flags `--amb-strand --sort-by-len --result 3` (the TR sequence in the reference is in red):

![](figures/109.png)
<!-- ![](figures/59.png) -->
<!-- ![](figures/60.png) -->




<!--
## Characterizing reads with clipped alignments

In the ST001 230x PacBio sample, 3.1% of all reads (i.e. 1'331'049 of 42'912'521 total) have at least one clipped alignment, defined as an alignment with either a 100bp left or a 100bp right soft clip and that involves a canonical chromosome (see [ClippedAlignments.java](https://github.com/fabio-cunial/smaht_experiments/blob/main/scripts/ClippedAlignments.java)). Here are some reads with clipped alignments, where colors corresponds to chromosomes in forward (positive) and reverse-complement (negative) orientation, and where:

* 0 = no alignment
* 23 = chrX
* 24 = chrY
* 25 = chrM
* 30 = non-canonical

![](figures/66.png)

Many reads map to just one chromosome and then continue with sequence that is not mapped to any region of GRCh38, and many reads map to at least one non-canonical GRCh38 chromosome:

![](figures/67.png)

The 230x PacBio BAM was created by merging the following 6 files provided in [the benchmarking section of the data portal](https://data.smaht.org/data/benchmarking/donor-st001#liver) and generated by different sequencing centers (an asterisk identifies Fiber-seq samples):

```
* ST001-1A-003S3-M22-B004-uwsc-SMAFIAK7UKIU-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
  ST001-1A-XX-M22-B001-bcm-SMAFI1DZ62NM-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
  ST001-1A-XX-M22-B001-bcm-SMAFI9VER4C3-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
  ST001-1A-XX-M22-B001-bcm-SMAFIYT6BC6F-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
  ST001-1A-XX-M22-B001-broad-SMAFIFXDVRTL-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
  ST001-1A-XX-M22-B001-washu-SMAFIQB3OD1L-pbmm2_1.13.0_GRCh38.aligned.sorted.bam
```

By a superficial IGV inspection, it seems that every one of these files contains reads with clipped alignments, e.g.:

![](figures/68.png)

By running the previous analysis on each BAM separately, we get the following counts:

```
Total reads: 4226911   Reads with clipped alignments: 129911 (3.07%)
Total reads: 4899244   Reads with clipped alignments: 129582 (2.64%)
Total reads: 4734888   Reads with clipped alignments: 124123 (2.62%)
Total reads: 5099541   Reads with clipped alignments: 130557 (2.56%)
Total reads: 12125667  Reads with clipped alignments: 346787 (2.86%)
Total reads: 11826270  Reads with clipped alignments: 470089 (3.97%)
```

From the SMaHT benchmarking manuscript, it seems that all sequencing centers got a similar mix of cells: 

_"A single lab prepared homogenized tissues and distributed samples to each Genome Characterization Center. Tissue homogenization ensured that the same genomic content was utilized across the centers to check for consistency and robustness."_ 

_"For each tissue, five 10 cm x 1 cm x 1 cm samples from defined anatomical locations in each organ were rinsed with saline, minced into several small pieces, placed in individual 50 ml conical tubes and snap frozen."_ 

_''To minimize regional variability and ensure consistency across sequencing centers, frozen tissues were pulverized into a powder using a pre-chilled mortar and pestle with liquid nitrogen. The homogenized material was pooled, mixed for uniformity, and aliquoted."_

Presence in other tissues/samples (PacBio only, by superficial IGV inspection only):

```
ST001	liver	Y
ST001	lung	Y
ST002	colon	Y
ST002	lung	Y
ST003	brain	Y
ST004	brain	Y
```

These seem to be even more frequent in ONT, but we have to make sure ONT was aligned to the same reference (see error message from header in CRAMs).
-->








# Genes associated with liver function

**Note: some interesting regions in or around these genes might be missing from what follows. Go over IGV in detail once again.**

We inspect each gene with IGV and we focus on regions that contain clusters of nearby SVs, clusters of clipped alignments, or clusters of SNPs with different patterns. In every IGV screenshot, the top track shows SMHT001 (cirrotic) and the bottom track shows ST001 (healthy). We perform two analyses:

* Do clipped alignments support different variants from spanning alignments? We try to answer this question by extracting all alignments in a window, from both ST001 and SMHT001, with the "Export alignments" feature of IGV. We convert them to FASTA with `samtools fasta -@ 16 file.sam` and we build a POA graph using abPOA in local alignment mode with flags `-m 1 --amb-strand --sort-by-len --result 3` (global alignment creates much more complex graphs). We show in blue POA nodes that are traversed by >=3 reads and in red nodes that belong to the reference sequence path.
* Do spanning reads in the ST001 230x PacBio BAM support >2 haplotypes? We try to answer this question by extracting just the spanning portion of each read, with the hapestry command that follows, and then by loading the FASTA into Jalview and using MAFFT (preset `FFT-NS-1`, speed oriented).

```
extract_reads_from_windows --output_dir ./reads_with_q/ \
	--bam_csv ${SAMPLES_LIST} \
	--windows ${FLANKED_WINDOWS_BED} \
	--bam_not_hardclipped \
	--require_spanning \
	--flank_length 0 \
	--fetch_max_length 100000 \
	--tags NM \
	--n_threads 1 \
	--force_forward
```

Often several clipped alignments seem to have mismatching bases aligned with a hom or het INS, and the patterns of mismatch seem to be supported by more than one read. Of course these may just be the beginning/end of the same INS. Often variation occurs inside TRs.




## From Ng et al. 2021

Ng, Stanley WK, et al. "[Convergent somatic mutations in metabolism genes in chronic liver disease.](https://www.nature.com/articles/s41586-021-03974-6)" Nature 598.7881 (2021): 473-478.

**Some notes on the paper.** Interestingly every one of their short-read samples is just 31x, but they do take multiple microdissections from the same liver, at controlled distances between microdissections (to probe different clones and prove convergent mutation). It's also interesting that they find structural variants at or near AVCR2A, GPAM, FOXO1, and that they estimate telomere lengths (another analysis for which long reads may be superior).

See directory [figures/ng_et_al_2021](https://github.com/fabio-cunial/smaht_experiments/tree/main/figures/ng_et_al_2021) for a full list of screenshots. In what follows we only show genes with potentially interesting features.


### ⭐️ ABCB11

[gnomAD SV >](https://gnomad.broadinstitute.org/region/2-168968743-168981935?dataset=gnomad_sv_r4)

![](figures/ng_et_al_2021/ABCB11.png)

<!-- #### Global alignment
![](figures/81.png)
The reads seem to support the hom INS, but they also create a more complex topology before/after it.
-->

The POA graph seems to support the presence of multiple INS:

![](figures/91.png)

By selecting some arbitrary nodes in each branch of the graph, and by enumerating all paths that use any of these nodes, there seem to be one main haplotype and a second one with smaller read support:

```
nReads   path
   106   28157+,28384+,28808+,28877+,28969+,29184+,29329+,29495+,29687+,
     6   29687-,29329-,29184-,28877-,28808-,28157-, 
```

![](figures/107.png)

However the MSA of all spanning reads seems to show two frequent haplotypes:

![](figures/94.png)



---

### B9D2

[gnomAD SV >](https://gnomad.broadinstitute.org/region/19-41351836-41379435?dataset=gnomad_sv_r4)

![](figures/ng_et_al_2021/B9D2.png)

<!-- #### Global alignment
![](figures/82.png)
The reads seem to support the hom INS, but they also create a more complex topology, including completely divergent paths.
-->

There seems to be just one INS:

![](figures/93.png)








## From the Goodman Lab

Genes with clear associations with liver disease, e.g. genes that are either known genetic risk factors for fatty liver disease, or protect from alcohol-related liver disease. The [Goodman lab](https://goodmanlab.mgh.harvard.edu/) has interests, tools and techniques to study GCKR, ADH1B, and MLXIPL, so any novel biology related to those three genes would have the lowest activation energy for mechanistic studies.

```
PNPLA3, TM6SF2, APOE, GCKR, TRIB1, GPAM, MARC1, MTTP, ADH1B, TOR1B, TMC4/MBOAT7, COBLL1, SREBF1, INSR, FTO, PNPLA2, MTARC1, MLXIPL, ADH1C, HFE, ATP7B, FRZB, IL18RAP, FLT3, GDF15, HGFAC, FSTL3, INHBA, INHBB
```

See directory [figures/goodman](https://github.com/fabio-cunial/smaht_experiments/tree/main/figures/goodman) for a full list of screenshots. In what follows we only show genes with potentially interesting features.
<!--
### GCKR
The gene seems to have no SV but >=2 methylation patterns:
![](figures/61.png)
![](figures/62.png)
------------------- OLD, delete what follows --------------------------
We extract all alignments (spanning and non-spanning) with the "Export alignments" feature of IGV. Then, we build a POA graph using abPOA with flags `-m 0 --amb-strand --sort-by-len --result 3` (the gene sequence in the reference is in red):

![](figures/63.png)
![](figures/64.png)
![](figures/65.png)
-->




### ⭐️ MTTP

[gnomAD SV >](https://gnomad.broadinstitute.org/region/4-99572824-99625997?dataset=gnomad_sv_r4)

![](figures/goodman/MTTP.png)
![](figures/103.png)

<!--#### Global alignment
![](figures/74.png)
The reads seem to support a discontinuity in the middle of the reference window, where the hom INS coexists with divergent paths.
-->

There seem to be multiple nearby INS, which may be the same ones described by the spanning alignments in the BAM:

![](figures/85.png)

By selecting some arbitrary nodes in each branch of the graph, and by enumerating all paths that use any of these nodes, there seem to be at least 3 well-supported haplotypes, which combine the 5 INS in different ways:

```
nReads   path
    66   16861+,17093+,17223+,17289+,17491+,17515+,17619+,17725+,17806+,
    70   16861+,17093+,17223+,17289+,17491+,17515+,17725+,17806+,
     4   17806-,77420-,17725-,17619-,17515-,17491-,17289-,77357-,17223-,17093-,16861-,
```

![](figures/102.png)

MSA of all spanning reads:

![](figures/98.png)
![](figures/99.png)

TRGT detects two haplotypes:
```
chr4	99588686	.	

AATATATATATACACACATATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATGTTCATATATATATTCATATATATATATATTTTTTTTTT	

AATATATATATACACACATATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATGTTCATATATATTCATATATATATGTTCATATATATTCATATATATGTTCATATATATATTCATATATATGTTCATATATATTCATATATATGTTCATATATATATTCATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATATATATTTTTTTT,

AATATATATATACACACATATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATATGTTCATATATATACACACATATATATGTTCATATATATACACACATATATATGTTCATATATATATTCATATATATATGTTCATGTATATTCATATATATATGTTCATATATATATTCATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATATGTTCATATATATATTCATATATATGTTCATATATATATTCATATATATATATATTTTTTTTTT
```
And 250 reads span them (coverage in the region: 300), but with no consistent CIGAR pattern.

Sniffles reports one 64bp INS:
```
chr4	99588927

A

AATGTTCATATATATATTCATATATATGTTCATATATATTCATATATATGTTCATATATATATTC
```

### ⭐️ TMC4

[gnomAD SV >](https://gnomad.broadinstitute.org/region/19-54168466-54169263?dataset=gnomad_sv_r4)

![](figures/goodman/TMC4.png)

There seem to be more than two haplotypes:

![](figures/71.png)

TRGT detects two haplotypes:
```
chr19	54168770	.	

CCTTTCTTTCTTTCTTTTCTTTTCTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTCTTTCCTTTCTTTCTTCTTTCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCTTCTTTCTCTCTCTCTCTCTCTCTCTCTATATATATATATATATTTTTTTTTTTTTCTTTTCTTTTCTTTTCTTTTTTTTTTTT	

CCTTTCTTTCTTTCTTTTCTTTTCTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTCTTTCCTTTCTTTCTTCTTTCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCTTCTTTCTCTCTCTCTCTCTCTCTCTCTATATATATATATATTTTTTTTTTTTTTCTTTTCTTTTCTTTTCTTTTTTTTTTTT,

CCTTTCTTTCTTTCTTTTCTTTCTTTTCTTTCTTTCTTTCTTTTCTTTCTTTCTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTCTTTCTTTCTTTCTTTCTTTCTTCCTTTCTTTCTTTTCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTTCTTTCTCTCTCTCTCTCTCTCTATATATATATATATATATTTTTCTTTTCTTTTCTTTTCTTTTTTTTTTTTT
```
And 123 reads span them (coverage in the region: 140). Aligning the reads to the two haplotypes shows some recurrent CIGAR edits, e.g.:

![](figures/120.png)

Sniffles: no support.


### ⭐️ MBOAT7

[gnomAD SV >](https://gnomad.broadinstitute.org/region/19-54181019-54182615?dataset=gnomad_sv_r4)

![](figures/goodman/MBOAT7.png)

There seem to be more than two haplotypes:

![](figures/70.png)

TRGT detects only one haplotype:
```
chr19	54181615	.	

CAGGGAGGGAGGGAGGGAGGGAGGGAAGGAGTGAAGAAGGGAAGAAAGAAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGAAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAGGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAAGGAAGGAAG	

CAGGGAGGGAGGGAGGGAGGGAGGGAAGGAGTGAAGAAGGGAAGAAAGAAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGAAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGAAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAAGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAGGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAGGGAAGGAGGGAAGGAAGGAGGGAGGGAAGGAAGGAGGGAAGGAAGGAAG
```
And 87 reads span it (coverage in the region: 100). Aligning the reads to the REF and ALT does not show recurrent CIGAR patterns:

![](figures/121.png)

Sniffles: no support.


---

### PNPLA3

[gnomAD SV >](https://gnomad.broadinstitute.org/region/22-43921805-43949582?dataset=gnomad_sv_r4)

This variation is unlikely to be somatic, but it is interesting because it affects the cirrotic sample but not the healthy sample. 

![](figures/goodman/PNPLA3.png)

<!--#### Global alignment
![](figures/72.png)
The reads seem to support a discontinuity in the middle of the reference window, with completely divergent paths.
-->

The variant seems to be a 1.2kb LINE1 INS:

![](figures/83.png)





### ADH1C

[gnomAD SV >](https://gnomad.broadinstitute.org/region/4-99334497-99354746?dataset=gnomad_sv_r4)

![](figures/goodman/ADH1C.png)

<!--#### Global alignment
![](figures/75.png)
-->

It seems that all left clips continue on chr20 or on a consistent chrUn, and that all right clips are short and unmapped:

![](figures/76.png)
![](figures/86.png)

The simplest explanation might be that chr20 contains a new copy of that intronic segment, and in fact the intronic segment is an L1:

![](figures/97.png)




### INSR

[gnomAD SV >](https://gnomad.broadinstitute.org/region/19-7127850-7178939?dataset=gnomad_sv_r4)

![](figures/goodman/INSR.png)
![](figures/101.png)

<!--#### Global alignment
![](figures/73.png)
The reads seem to support a discontinuity in the middle of the reference window, where the hom INS coexists with divergent paths.
-->

There seem to be more than one INS, which may be the same ones described by the spanning alignments in the BAM:

![](figures/84.png)

By selecting some arbitrary nodes in each branch of the graph, and by enumerating all paths that use any of these nodes, there seem to be mainly two haplotypes, with a third hap. only supported by one read:

```
nReads   path
    45   10542-,10389-,10198-,10014-,9858-,9690-,
    45   10542-,64302-,10389-,64155-,63933-,10198-,10014-,9858-,9690-,
     1   10542-,10389-,10198-,9858-,9690-,
```

![](figures/100.png)

MSA of all spanning reads:

![](figures/95.png)
![](figures/96.png)

TRGT detects two haplotypes:
```
chr19	7152980	.	

TACACACACACACACACCCCACACACACACACACCACACACACACACCACACACACACACACCACACACCCCCCCACACACACACACCACACACATCACACAACACACCACACATACACACACCACAAACACACAACCACACACACACACCACACACCACACACACCACACAACACACCACACATACACACACCACAAACACACAACCACACACACACCACACACACCACACACCACACACACGCACCACACACACACCATACACGCACCACACACACACTACATGCACACCTCACACACCACACCACACACCGCACACACACCACACAAATCACACACACACCACACACACCCACGCCACACACCACAGACCACACACCACACACCCCACACACACCATACACGCACCACACACACACCACATACACACTTCACACACACCACACACCACACAAATCACACACACCACCACACCCACGCCACACACCACACACACACCACACACGCACCACACACACACACAGTGAAG	

TACACACACACACACACCCCCCACACACACACACCACACACACACACCACACACACCACACATACACACACCACACACCCCCCCCACACACACACCACACACATCACACAACACACCACACACCACACATACACACACCACAAACACACAACCACACACACACACCACACACACCACACACACCACACAACACACCACACATACACACACCACAAACACACAACCACACACACCACACACACCACACACCACACACACCACACAACACACCACACATACACACACCACAAACACAACCACACACACACACCACACACACCACACACACCACACACGCACCACACACACACCATACACGCACCACACACACACTACATGCACACCTCACACACCACACCACACACCGCACACACACCACACAAATCACACACACACCACACACACCCACGCCACACACCACAGACCACACACACCACACACCCCACACACACCACACACCACACACACCAAACACACCATACATGCACCACACACACACCACATGCACACCTCACACACCATACCACACACCACACACACCACACAAATCACACACATACTACACACACCCACGCCACACACCACACACACACCACACACGCCACACAACGGACCACACATACACCACACACCACAAACACAACCACACACACACTATACACACACACCACACACACCACACAACAGACCACACATACACACACCACACACCACAAACACACACCATACACACACACCACACATACCACACACACCACACACACAACCACACACACACACCACACACACAACACACAACAGACCACACATACACACACCACAAACACAACCACACACACCACACATACCACACACACACCACACATACCACACACACACCACACACACAACCACACACACCACACCATACACACAACCACACACACACCATACACGCACCACACACACACCACATACACACTTCACACACACCACACACCACACAAATCACACACACCACCACACCCACGCCACACACCACACACACACCACACACGCACCACACACACACACAGTGAAG,

TACACACACACACACCCCCCCCACACACACACACCACACACACACACCACACACACCACACATACACACACCACACACCCCCCCACACACACACCACACACATCACACAACACACCACACATACACACACCACAAACACACAACCACACACACACACCACACACCACACACACCACACAACACACCACACATACACACACCACAAACACACAACCACACACACACCACACACCACACACCACACACACCACACAACACACCACACATACACACACCACAAACACACAACCACACACACACACCACACACACCACACACACCACACACGCACCACACACACACCATACACGCACCACACACACACTACATGCACACCTCACACACCACACCACACACCGCACACACACCACACAAATCACACACACACCACACACACCCACGCCACACACCACAGACCACACACACCACACACCCCACACACACCACACACCACACACACCAAACACACCATACATGCACCACACACACACCACATGCACACCTCACACACCATACCACACACCACACACACCACACAAATCACACACATACTACACACACCCACGCCACACACCACACACCACACACGCCACACAACGGACCACACATACACCACACACCACAAACACAACCACACACACACTATACACACACACCACACACACCACACAACAGACCACACATACACACACCACACACCACACACACCATACACACACACCACACATACCACACACACCACACACACAACCACACACACACACCACACACACAACACACAACAGACCACACATACACACACCACAAACACAACCACACACACCACACCATACACACACACCATACACGCACCACACACACACCACATACACACTTCACACACACCACACACCACACAAATCACACACACCACCACACCCACGCCACACACCACACACACACCATACACGCACCACACACACACTACATGCACACCTCACACACCACACCACACACCGCACACACACCACACAAATCACACACACACCACACACACCCACGCCACACACCACAGACCACACACACCACACACCCCACACACACCACACACCACACACACCAAACACACCATACATGCACCACACAACCACATGCACACCTCACACACCATACCACACACCACACACACCACACAAATCACACACATACTACACACACCCACGCCACACACCACACACACACCACACACGCCACACAACGGACCACACATACACCACACACCACAAACACAACCACACACACACTATACACACACACCACACACACCACACAACAGACCACACATACACACACCACACACCACACACACCATACACACACACCACACATACCACACACACCACACACACAACCACACACACACACCACACACACAACACACAACAGACCACACATACACACACCACAAACACAACCACACACACCACACCATACACACACACCATACACGCACCACACACACACCACATACACACTTCACACACACCACACACCACACAAATCACACACACCACCACACCCACGCCACACACCACACACGCACCACACACACACACAGTGAAG
```
And 135 reads span them (coverage in the region: 160). Aligning the reads to the two haplotypes does not show recurrent CIGAR patterns.

Sniffles report one 77bp INS:
```
chr19	7153186

A

CCACACACACACCACACACCACACACCACACACACCACACAACACACCACACATACACACACCACAAACACACAACCA
```


### FLT3

[gnomAD SV >](https://gnomad.broadinstitute.org/region/13-27987701-28012619?dataset=gnomad_sv_r4)

![](figures/goodman/FLT3.png)
![](figures/104.png)

<!--#### Global alignment
![](figures/80.png)
There seem to be 3 hom INS as in the BAM. However, one INS is located just outside the region, and tracing some reads with clipped alignments in the graph gives some unexpected results WRT the BAMs, so I'm not confident about the output of abPOA in this window.
-->

There seem to be just 3 INS, which may be the same ones described by spanning alignments in the BAM:

![](figures/90.png)

By selecting some arbitrary nodes in each branch of the graph, and by enumerating all paths that use any of these nodes, there seems to be just one haplotype:

```
nReads   path
    74   1268-,730-,312-,158-,103-,71038-,70873-,70683-,70611-,
```

![](figures/106.png)

However, the MSA of all spanning reads seems to show two haplotypes:

![](figures/105.png)

TRGT detects two haplotypes:
```
chr13	27998700	.	

ATGTGTGTGTGTGTGTATATGTAACTATGTGTGTGCAGTTATATATATATAAATATATATATTTATATATAAATATATGTTATATATACATAAATATATGTTATATATATTTATATATATTTATATACATATAAATATATATATTTATATACATATATAAAAATATATATATATTTATATATAAATATATGTTATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATTTATATACATATAAATATATGTATATACATATAAATATATATTTATATACATATAAATATATATTTATATACATATAAATATATGTATATACATATATATATGTATATACATATAAATATATATGTATATACATATAAATATATATGTATATACATATAAATATATATAATATTTATATACATATATATATTTATATTTATATGTATATACATATATATTTATATTTATATGTATATACATATATATATTTATATTTATATGTATATACATATATATATGTATATAAATATATATATTTATATGTATATAATATACATACATATATATATAT	ATGTGTGTGTGTGTGTATATGTAACTATGTGTGTGCAGTTATATATATATAAATATATATATTTATATATAAATATATGTTATATATATAAATATATGTTATATATATAAATATATGTTATATATATAAATATATGTTATATATATAAATATATGTTATATATATAAATATATGTTATATATATTTATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATATAAATATATGTTATATATATTTATATGCATATAAATATATATTTATATGCATATAAATATATATTTATATGCATATAAATATATATTTATATACATATAAATATATATGTATATACATATAAATATATGTATATACATATAAATATATATGTATATACATATATATATGTATATACATATATATATATGTATATACATATAAATATATATGTATATACATATAAATATATGTATATACATATAAATATATATGTATATACATATAAATATATATGTATATACATATAAATATATATGTATGTATATTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTACATACATATAAATATATATTTACATACATATAAATATATACATTTATATACATATAAATATATACATATATACATATAAATATATATTTATATACATATAAATATATATATTTATATATATTTATATATATATAAATATATATATATTTGTATATATATATAAATATATATATATTTATATATATATAAATATATATATATTTATATATATATACATATAAATATATATTTATATACATATAAATATATAATTATATACATATAAATATATAATTATATACATATAAATATATATAATTATATACATATAAATATATATTTATATACATATAAATTTATATACATATAAATATAAATATTTATATTTATATACATATAAATATATATATTTATATTTATATACATATATATTTATATTTATATGTATATACATATATATATTTATATTTATATGTATATACATATATATATTTATATTTATGTATATACATATATATATTTATATTTATATGTATATAAATATATATATTTATATTTATATAATATACATACATATATATATAT,

ATGTGTGTGTGTGTGTATATGTAACTATGTGTGTGCAGTTATATATATATAAATATATATATTTATATATAAATATATGTTATATATATAAATATATGTTATATATATAAATATATGTTATATATATAAATATATATATTTATATATAAATATATGTTATATATATAAATATATGTTATATATATTTATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATATAAATATATGTTATATATATTTATATGCATATAAATATATATTTATATGCATATAAATATATATTTATATGCATATAAATATATATTTATATACATATAAATATATATGTATATACATATAAATATATGTATATACATATAAATATATATGTATATACATATAAATATATATGTATATACATATATATATATGTATATACATATAAATATATATGTATATACATATATATATGTATATACATATATATATATGTATATACATATAAATATATATGTATATACATATAAATATATGTATATACATATAAATATATATGTATATACATATAAATATATATGTATATACATATAAATATATATGTATGTATATTATATACATATAAATATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTATATACATATAAATATATATATTTATATACATATAAATATATATTTACATACATATAAATATATATTTACATACATATAAATATATACATTTATATACATATAAATATATACATATATACATATAAATATATATTTATATACATATAAATATATATATTTATATATATTTATATATATATAAATATATATATATTTGTATATATATATAAATATATATATATTTATATATATATAAATATATATATATTTATATATATATACATATAAATATATATTTATATACATATAAATATATAATTATATACATATAAATATATAATTATATACATATAAATATATATAATTATATACATATAAATATATATTTATATACATATAAATTTATATACATATAAATATAAATATTTATATTTATATACATATAAATATATATATTTATATTTATATACATATATATTTATATTTATATGTATATACATATATATATTTATATTTATATGTATATACATATATATATTTATATTTATGTATATACATATATATATTTATATTTATATGTATATAAATATATATATTTATATTTATATAATATACATACATATATATATAT
```
And 138 reads span them (coverage in the region: 150). It's unclear if aligning the reads to the two haplotypes shows recurrent CIGAR patterns, e.g.:

![](figures/122.png)

Sniffles: no support.


### ATP7B

![](figures/goodman/ATP7B.png)

<!--#### Global alignment
![](figures/77.png)
It seems that the clipped alignments support the hom INS. However, abPOA does not position the INS at the center of the reference window, and tracing some reads with clipped alignments in the graph gives some unexpected results WRT the BAMs, so I'm not confident about the output of abPOA in this window.
-->

There seems to be just one INS, as described by spanning alignments in the BAM:

![](figures/87.png)




### IL18RAP

![](figures/goodman/IL18RAP.png)

<!--#### Global alignment
![](figures/78.png)
It seems that the clipped alignments support the hom INS. However, tracing some reads with clipped alignments, or some spanning reads with the INS, in the graph, gives some unexpected results WRT the BAMs, so I'm not confident about the output of abPOA in this window.
-->

There seems to be just one INS, as described by spanning alignments in the BAM:

![](figures/88.png)




### SREBF1

![](figures/goodman/SREBF1.png)

<!--#### Global alignment
![](figures/79.png)
Rather than an INS, there seem to be two alternative sequences.
-->

There seems to be just one INS, as described by spanning alignments in the BAM:

![](figures/89.png)




### HGFAC

![](figures/goodman/HGFAC.png)

There seem to be just two haplotypes:

![](figures/69.png)

TRGT detects two haplotypes:
```
chr4	3440598	.	

TCCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCCCCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCC	

TCCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCC,
TCCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCCCCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCATCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACTC

.	.	

TRID=4-3440598-3440776-CCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCC;END=3440814;MOTIFS=CCCAGCACCCACCATGCCTCCACATTCACCGTGCCCCCAGGGCCCACCATGACCTCCCCCCCAACCCCCAACCCCCCATTCACCGTGCC;STRUC=<VC172810>

GT:AL:ALLR:SD:MC:MS:AP:AM	

1/2:74,323:71-75,317-325:73,59:1,5:0(0-74),0(0-323):0.831461,0.723596:0.65,0.63
```
And 132 reads span them (coverage in the region: 140). Aligning the reads to the two haplotypes seems to show some recurrent CIGAR operations, e.g.:

![](figures/123.png)

Sniffles: no support.


## From Brzozowska et al. 2025

Brzozowska, N. et al. "[Selection for somatic escape variants in SERPINA1 in the liver of patients with alpha-1 antitrypsin deficiency.](https://www.nature.com/articles/s41588-025-02125-1)" Nature Genetics (2025): 1-9.

### SERPINA1

[gnomAD SV >](https://gnomad.broadinstitute.org/region/14-94361139-94366969?dataset=gnomad_sv_r4)

![](figures/111.png)
![](figures/112.png)

There seems to be just the 2kb hom INS described by the BAM:

![](figures/113.png)

It might be another copy of SERPINA2:

![](figures/114.png)


---


# SAVANA

We tried to run [SAVANA](https://www.nature.com/articles/s41592-025-02708-0) [1.3.6](https://github.com/cortes-ciriano-lab/savana) on the 230x PacBio ST001 liver. Specifically, we ran `savana to --pb --tumour ~{aligned_bam} --ref ~{reference_fa} --g1000_vcf 1000g_hg38 --contigs /savana/example/contigs.chr.hg38.txt`. Remarks:
* SAVANA seems to be designed for tumor-normal input, but it has a tumor-ony mode (`to`) and we used that in the experiment.
* It seems to be designed for ONT, but it has a PacBio mode (`--pb`) and we used that in the experiment.

Unfortunately it crashed after 6.5h as follows:
```
[...]
INFO:root:Breakpoint 25663 rejected
INFO:root:Testing validity of 25716 in interval from 7393 to 25822 yields statistic 3.4292046763790474
INFO:root:Breakpoint 25716 accepted
Traceback (most recent call last):
  File "/opt/conda/bin/savana", line 8, in <module>
    sys.exit(main())
             ^^^^^^
  File "/opt/conda/lib/python3.11/site-packages/savana/savana.py", line 762, in main
    args.func(args)
  File "/opt/conda/lib/python3.11/site-packages/savana/savana.py", line 367, in savana_tumour_only
    savana_cna(args, True)
  File "/opt/conda/lib/python3.11/site-packages/savana/savana.py", line 277, in savana_cna
    fit_absolute.fit_absolute_cn(outdir, log2r_cn_path, allele_counts_bed_path, args.sample,
  File "/opt/conda/lib/python3.11/site-packages/savana/fit_absolute.py", line 167, in fit_absolute_cn
    fits = cnfitter.estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=weights, distance_function=distance_function)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/conda/lib/python3.11/site-packages/savana/cn_functions.py", line 233, in estimate_grid_distances
    d = acn_distance(relative_CN, cur_pur, cur_ploi, weights=weights, distance_function=distance_function)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/conda/lib/python3.11/site-packages/savana/cn_functions.py", line 181, in acn_distance
    acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/conda/lib/python3.11/site-packages/savana/cn_functions.py", line 181, in <listcomp>
    acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/conda/lib/python3.11/site-packages/savana/cn_functions.py", line 173, in relative_to_absolute_CN
    acn = ploidy + (relative_CN - 1)*(ploidy+(2/purity)-2) 
                                              ~^~~~~~~
ZeroDivisionError: division by zero
```

On a VM with 32 cores and 64GB of RAM, it used 28 cores and 20GB of RAM, and it cost $10.


---

# Telomere length analysis

See this [concept sheet](https://docs.google.com/document/d/1T1aTPMkV_C-9Z-zc1X7TYbYaD1YWzG5kRMGZUFfhb5o/edit?tab=t.0). This person has already developed a tool and a paper, and is joining Tim's lab soon.

Cells under stress might replicate more often. TERT is the telomerase gene, and mutations in its promoter might affect telomere length. E.g. mutations could lead to liver cancer if they allow to escape checks on telomere length (Natalia).




# Cell type deconvolution

Separate reads from mature hepatocytes vs progenitors. Using methylation or Fiber-seq. May depend on how good the annotation of the cell types we want is. However we could do it reference-free, by just clustering the methylation patterns of the reads.




# Combining tissues from the same donor to detect embryonic variants

All 60x PacBio samples should be available in the data portal (in the other section than benchmarking). But the portal is down now.



# Donors SV paper

Joining the SV donors paper makes sense, esp. for doing a deeper dive on some genes or types of variation.



# Long-read RNA

Might be useful also for detecting global splicing patterns (i.e. systematic differences in which transcripts are expressed), induced by a global splicing cofactor gene (Natalia).



---

# References to read

https://pmc.ncbi.nlm.nih.gov/articles/PMC12339610/

https://trtools.readthedocs.io/en/stable/source/prancSTR.html

https://academic.oup.com/bioinformatics/article/40/8/btae485/7723996

TR and liver: https://research.edgehill.ac.uk/en/publications/length-of-variable-numbers-of-tandem-repeats-in-the-carboxyl-este#:~:text=Length%20of%20Variable%20Numbers%20of%20Tandem%20Repeats,M%C3%B6ssner%2C%20Claudia%20Ruffert%2C%20Mario%20Krehan%2C%20Christian%20Zapf

Which liver genes have TRs in exons?

https://pubmed.ncbi.nlm.nih.gov/40004542/

https://www.mdpi.com/2073-4425/16/2/213


RNA aberrant splicing:
https://github.com/gagneurlab/fraser?tab=readme-ov-file

## From Natalia

https://www.nature.com/articles/s41586-019-1670-9

https://www.nature.com/articles/s41586-021-03974-6

https://pubmed.ncbi.nlm.nih.gov/31597092/

https://www.nature.com/articles/s41467-024-55325-4



# More on hepatocytes

Hepatocites become tetraploid at some point in their lifetime (Natalia).

Different types of hepatocyte:
* Zone 1 (peri-portal vein)
* Zone2
* Zone3 (peri-central vein)
There is an oxygen gradient, since blood comes in with a lot of oxygen, and the hepatocytes in the three zones might have slightly diffeent metabolic function. However, they do not seem to separate in RNA umaps.

If there is rewiring of some metabolic pathways due to somatic mutations, in which zone is it happening? Maybe we could use methylation to assign variants to zones.

Other cell types in the liver:
* Cholangiocytes
* Kupfer cells.
