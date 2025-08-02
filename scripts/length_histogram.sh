#!/bin/bash
#
INPUT_VCF="../vcfs/ST001_230_sniffles.vcf.gz"
REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
LENGTHS="100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000"

# ST001: length histogram at full coverage.
rm -f ST001_230_lengths.csv
PREVIOUS_LENGTH="0"
for LENGTH in ${LENGTHS}; do
    TOTAL_CALLS=$(bcftools view --no-header --regions ${REGIONS} ${INPUT_VCF} | wc -l)
    INS_CALLS=$(bcftools filter --include "SVTYPE=INS && ((SVLEN>${PREVIOUS_LENGTH} && SVLEN<=${LENGTH}) || (SVLEN<-${PREVIOUS_LENGTH} && SVLEN>=-${LENGTH}))" --regions ${REGIONS} ${INPUT_VCF} | wc -l)
    DEL_CALLS=$(bcftools filter --include "SVTYPE=DEL && ((SVLEN>${PREVIOUS_LENGTH} && SVLEN<=${LENGTH}) || (SVLEN<-${PREVIOUS_LENGTH} && SVLEN>=-${LENGTH}))" --regions ${REGIONS} ${INPUT_VCF} | wc -l)
    echo "${LENGTH},${TOTAL_CALLS},${INS_CALLS},${DEL_CALLS}" >> ST001_230_lengths.csv
    PREVIOUS_LENGTH=${LENGTH}
done

# SMHT001: number of calls.
TOTAL_CALLS=$(bcftools view --no-header --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l)
INS_CALLS=$(bcftools filter --include "SVTYPE=INS" --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l)
DEL_CALLS=$(bcftools filter --include "SVTYPE=DEL" --regions ${REGIONS} ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l)
INSIDE_TR=$(bcftools view --no-header --regions-file ${TR_BED} --regions-overlap pos ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l)
OUTSIDE_TR=$(bcftools view --no-header --regions-file nontr.bed --regions-overlap pos ${INPUT_DIR}/SMHT001_010_sniffles.vcf.gz | wc -l)
echo "${COVERAGE},${TOTAL_CALLS},${INS_CALLS},${DEL_CALLS},${INSIDE_TR},${OUTSIDE_TR}" > SMHT001_counts.csv

# ST001: number of calls.
rm -f ST001_counts.csv
for COVERAGE in $(seq 10 10 ${MAX_COVERAGE}); do
    TOTAL_CALLS=$(bcftools view --no-header --regions ${REGIONS} ${INPUT_DIR}/ST001_${COVERAGE}_sniffles.vcf.gz | wc -l)
    INS_CALLS=$(bcftools filter --include "SVTYPE=INS" --regions ${REGIONS} ${INPUT_DIR}/ST001_${COVERAGE}_sniffles.vcf.gz | wc -l)
    DEL_CALLS=$(bcftools filter --include "SVTYPE=DEL" --regions ${REGIONS} ${INPUT_DIR}/ST001_${COVERAGE}_sniffles.vcf.gz | wc -l)
    INSIDE_TR=$(bcftools view --no-header --regions-file ${TR_BED} --regions-overlap pos ${INPUT_DIR}/ST001_${COVERAGE}_sniffles.vcf.gz | wc -l)
    OUTSIDE_TR=$(bcftools view --no-header --regions-file nontr.bed --regions-overlap pos ${INPUT_DIR}/ST001_${COVERAGE}_sniffles.vcf.gz | wc -l)
    echo "${COVERAGE},${TOTAL_CALLS},${INS_CALLS},${DEL_CALLS},${INSIDE_TR},${OUTSIDE_TR}" >> ST001_counts.csv
done
