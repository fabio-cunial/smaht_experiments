version 1.0


#
workflow Manta {
    input {
        String sample_id
        File input_alignments
        File input_alignments_index
        File baseline_alignments
        File baseline_alignments_index
        
        File reference_fa
        File reference_fai
        
        Int n_cores
        Int ram_gb
    }
    parameter_meta {
        input_alignments: "BAM or CRAM. Fed to the `--tumorBam` flag."
        baseline_alignments: "BAM or CRAM. Fed to the `--normalBam` flag."
    }

    call MantaImpl {
        input:
            sample_id = sample_id,
            input_alignments = input_alignments,
            input_alignments_index = input_alignments_index,
            baseline_alignments = baseline_alignments,
            baseline_alignments_index = baseline_alignments_index,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            n_cores = n_cores,
            ram_gb = ram_gb
    }

    output {
         File output_vcf_gz = MantaImpl.output_vcf_gz
         File output_tbi = MantaImpl.output_tbi
    }
}


# Remark: Manta is run with default parameters. No crucial parameter seems to
# be exposed on the command line.
#
# Performance on a machine with 16 cores and 16 GB of RAM:
#
# COVERAGE  CPU%    TIME    RAM
# 120x
#
task MantaImpl {
    input {
        String sample_id
        File input_alignments
        File input_alignments_index
        File baseline_alignments
        File baseline_alignments_index
        
        File reference_fa
        File reference_fai
        
        Int n_cores
        Int ram_gb
    }
    parameter_meta {
        input_alignments: "BAM or CRAM. Fed to the `--tumorBam` flag."
        baseline_alignments: "BAM or CRAM. Fed to the `--normalBam` flag."
    }
    
    Int disk_size_gb = ceil(size(input_alignments, "GB")) + ceil(size(baseline_alignments, "GB")) + ceil(size(reference_fa, "GB")) + 200
    String docker_dir = "/smaht_experiments"

    command <<<
        set -euxo pipefail
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        df -h
        
        mkdir ./input/
        mv ~{input_alignments} ./input/
        mv ~{input_alignments_index} ./input/
        mkdir ./baseline/
        mv ~{baseline_alignments} ./baseline/
        mv ~{baseline_alignments_index} ./baseline/
        ls -laht
        tree
        INPUT_FILE=$( find ./input/ -type f -maxdepth 1 -name "*.bam" -o -name "*.cram" )
        BASELINE_FILE=$( find ./baseline/ -type f -maxdepth 1 -name "*.bam" -o -name "*.cram" )
        ${TIME_COMMAND} python2 ~{docker_dir}/manta/bin/configManta.py --tumorBam ${INPUT_FILE} --normalBam ${BASELINE_FILE} --referenceFasta ~{reference_fa} --runDir ./manta
        ${TIME_COMMAND} python2 ./manta/runWorkflow.py -j ${N_THREADS}
        ls -laht
        tree
        ${TIME_COMMAND} bcftools sort -m $(( ~{ram_gb} - 2 ))G --output-type z --output ~{sample_id}_manta.vcf.gz ./manta/results/variants/somaticSV.vcf.gz
        ${TIME_COMMAND} tabix -f ~{sample_id}_manta.vcf.gz
    >>>

    output {
        File output_vcf_gz = sample_id + "_manta.vcf.gz"
        File output_tbi = sample_id + "_manta.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
