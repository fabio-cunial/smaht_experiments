#!/bin/bash
#
INPUT_DIR="../vcfs"
MAX_COVERAGE="230"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
BCFTOOLS_COMMAND="/Users/fcunial/git/bcftools-1.22/bcftools"

set -euxo


for CHR in ${CHROMOSOMES}; do
    # SMHT001
    TOTAL_CALLS=$(${BCFTOOLS_COMMAND} view --no-header ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz ${CHR} | wc -l | xargs)
    echo "10,${TOTAL_CALLS}" > SMHT001_counts_${CHR}.csv
    # ST001
    rm -f ST001_counts_${CHR}.csv
    for COVERAGE in $(seq 10 10 ${MAX_COVERAGE}); do
        COV=$(printf "%03d" ${COVERAGE})
        TOTAL_CALLS=$(${BCFTOOLS_COMMAND} view --no-header ${INPUT_DIR}/ST001_${COV}_sniffles.vcf.gz ${CHR} | wc -l | xargs)
        echo "${COVERAGE},${TOTAL_CALLS}" >> ST001_counts_${CHR}.csv
    done
done
