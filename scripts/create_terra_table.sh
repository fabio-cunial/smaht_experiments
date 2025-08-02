#!/bin/bash
#
TABLE_NAME="liver_experiment_downsampling"
SAMPLE_ID="ST001"
TISSUE_ID="liver_1a"
PLATFORM="PacBio Revio"
DOWNSAMPLED_PREFIX="gs://fc-82113812-f2f7-416b-9c15-baeb9e98781a/bams/${SAMPLE_ID}/${TISSUE_ID}/coverages_pacbio/merged_pacbio_"

echo -e "entity:${TABLE_NAME}_id\tliver_state\tdownsampled_bai\tdownsampled_bam\tcoverage_approx\tsource_platform"
echo -e "SMHT001_010\tliver failure, alcohol abuse\tgs://fc-3756f59c-37a9-431b-9959-6601943d887a/results/PBFlowcell/m84037_250123_221103_s3/reads/ccs/aligned/m84037_250123_221103_s3.hifi_reads.bc2012.bam.bai\tgs://fc-3756f59c-37a9-431b-9959-6601943d887a/results/PBFlowcell/m84037_250123_221103_s3/reads/ccs/aligned/m84037_250123_221103_s3.hifi_reads.bc2012.bam\t10.616978806227257\t${PLATFORM}"
for COVERAGE in $(seq 10 10 550); do
    echo -e "${SAMPLE_ID}_$(printf "%03d" ${COVERAGE})\thealthy\t${DOWNSAMPLED_PREFIX}${COVERAGE}.bam.bai\t${DOWNSAMPLED_PREFIX}${COVERAGE}.bam\t${COVERAGE}\t${PLATFORM}"
done
