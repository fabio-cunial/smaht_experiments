#!/bin/bash
#
INPUT_DIR="../vcfs"
MAX_COVERAGE="230"
REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
TR_BED="/Users/fcunial/Downloads/human_GRCh38_no_alt_analysis_set.trf.sorted.bed"
REFERENCE_FAI="/Users/fcunial/Downloads/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
BCFTOOLS_COMMAND="/Users/fcunial/git/bcftools-1.22/bcftools"

set -euxo


bedtools complement -i ${TR_BED} -g ${REFERENCE_FAI} | grep -v 'random' | grep -v 'chrUn' | grep -v 'chrEBV' > nontr.bed

# SMHT001
TOTAL_CALLS=$(${BCFTOOLS_COMMAND} view --no-header --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l | xargs)
INS_CALLS=$(${BCFTOOLS_COMMAND} filter --include 'SVTYPE="INS"' --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | grep -v ^# | wc -l | xargs)
DEL_CALLS=$(${BCFTOOLS_COMMAND} filter --include 'SVTYPE="DEL"' --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | grep -v ^# | wc -l | xargs)
BND_CALLS=$(${BCFTOOLS_COMMAND} filter --include 'SVTYPE="BND"' --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | grep -v ^# | wc -l | xargs)
INSIDE_TR=$(${BCFTOOLS_COMMAND} view --no-header --regions-file ${TR_BED} --regions-overlap pos ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l | xargs)
OUTSIDE_TR=$(${BCFTOOLS_COMMAND} view --no-header --regions-file nontr.bed --regions-overlap pos ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l | xargs)
echo "10,${TOTAL_CALLS},${INS_CALLS},${DEL_CALLS},${BND_CALLS},${INSIDE_TR},${OUTSIDE_TR}" > SMHT001_counts.csv

# ST001
rm -f ST001_counts.csv
for COVERAGE in $(seq 10 10 ${MAX_COVERAGE}); do
    COV=$(printf "%03d" ${COVERAGE})
    TOTAL_CALLS=$(${BCFTOOLS_COMMAND} view --no-header --regions ${REGIONS} ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz | wc -l | xargs)
    INS_CALLS=$(${BCFTOOLS_COMMAND} filter --include 'SVTYPE="INS"' --regions ${REGIONS} ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz | grep -v ^# | wc -l | xargs)
    DEL_CALLS=$(${BCFTOOLS_COMMAND} filter --include 'SVTYPE="DEL"' --regions ${REGIONS} ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz | grep -v ^# | wc -l | xargs)
    BND_CALLS=$(${BCFTOOLS_COMMAND} filter --include 'SVTYPE="BND"' --regions ${REGIONS} ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz | grep -v ^# | wc -l | xargs)
    INSIDE_TR=$(${BCFTOOLS_COMMAND} view --no-header --regions-file ${TR_BED} --regions-overlap pos ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz | wc -l | xargs)
    OUTSIDE_TR=$(${BCFTOOLS_COMMAND} view --no-header --regions-file nontr.bed --regions-overlap pos ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz | wc -l | xargs)
    echo "${COVERAGE},${TOTAL_CALLS},${INS_CALLS},${DEL_CALLS},${BND_CALLS},${INSIDE_TR},${OUTSIDE_TR}" >> ST001_counts.csv
done
