#!/bin/bash
#
INPUT_VCF="../vcfs/ST001_230_sniffles.vcf.gz"
REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
LENGTHS="100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000"
BCFTOOLS_COMMAND="/Users/fcunial/git/bcftools-1.22/bcftools"

set -euxo


# ST001: length histogram at max coverage.
rm -f ST001_230_lengths.csv
PREVIOUS_LENGTH="0"
for LENGTH in ${LENGTHS}; do
    TOTAL_CALLS=$(${BCFTOOLS_COMMAND} filter --include "((SVLEN>${PREVIOUS_LENGTH} && SVLEN<=${LENGTH}) || (SVLEN<-${PREVIOUS_LENGTH} && SVLEN>=-${LENGTH}))" --regions ${REGIONS} ${INPUT_VCF} | grep -v ^# | wc -l | xargs)
    INS_CALLS=$(${BCFTOOLS_COMMAND} filter --include "SVTYPE=\"INS\" && ((SVLEN>${PREVIOUS_LENGTH} && SVLEN<=${LENGTH}) || (SVLEN<-${PREVIOUS_LENGTH} && SVLEN>=-${LENGTH}))" --regions ${REGIONS} ${INPUT_VCF} | grep -v ^# | wc -l | xargs)
    DEL_CALLS=$(${BCFTOOLS_COMMAND} filter --include "SVTYPE=\"DEL\" && ((SVLEN>${PREVIOUS_LENGTH} && SVLEN<=${LENGTH}) || (SVLEN<-${PREVIOUS_LENGTH} && SVLEN>=-${LENGTH}))" --regions ${REGIONS} ${INPUT_VCF} | grep -v ^# | wc -l | xargs)
    echo "${LENGTH},${TOTAL_CALLS},${INS_CALLS},${DEL_CALLS}" >> ST001_fullcov_lengths.csv
    PREVIOUS_LENGTH=${LENGTH}
done

# ST001: AF histogram at max coverage.
${BCFTOOLS_COMMAND} query --format '%INFO/SUPPORT,%INFO/VAF\n' --regions ${REGIONS} ${INPUT_VCF} > ST001_fullcov_support_vaf.csv
