version 1.0


# Remark: this workflow does not output germline SVs.
#
workflow Sniffles {
    input {
        String sample_id
        File input_bam
        File input_bai
        Float coverage
        
        File reference_fa
        File reference_fai
        File tandems_bed
        
        Int min_sv_length = 50
        Int min_supporting_reads = 2
        
        Int n_cores
        Int ram_gb
    }
    parameter_meta {
        coverage: "Estimate of the coverage of the BAM."
        min_supporting_reads: "Minimum number of reads that support a call."
    }

    call SnifflesImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            coverage = coverage,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            tandems_bed = tandems_bed,
            min_sv_length = min_sv_length,
            min_supporting_reads = min_supporting_reads,
            n_cores = n_cores,
            ram_gb = ram_gb
    }

    output {
         File output_vcf_gz = SnifflesImpl.output_vcf_gz
         File output_tbi = SnifflesImpl.output_tbi
    }
}


# Performance on a machine with 16 cores and 32 GB of RAM:
#
# COVERAGE  CPU%    TIME    RAM
# 10x
#
task SnifflesImpl {
    input {
        String sample_id
        File input_bam
        File input_bai
        Float coverage
        
        File reference_fa
        File reference_fai
        File tandems_bed
        
        Int min_sv_length
        Int min_supporting_reads
        
        Int n_cores
        Int ram_gb
    }
    parameter_meta {
        coverage: "Estimate of the coverage of the BAM."
        min_supporting_reads: "Minimum number of reads that support a call."
    }
    
    Int disk_size_gb = ceil(size(input_bam, "GB")) + ceil(size(reference_fa, "GB")) + 200
    String docker_dir = "/smaht_experiments"

    command <<<
        set -euxo pipefail
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        df -h
        
        MOSAIC_AF_MIN=$(echo "scale=4; ~{min_supporting_reads}/~{coverage}" | bc)
        mkdir ./sniffles_tmp
        ${TIME_COMMAND} sniffles --threads ${N_THREADS} --tmp-dir ./sniffles_tmp \
            --mosaic \
            --mosaic-af-min ${MOSAIC_AF_MIN} \
            --output-rnames \
            --minsvlen ~{min_sv_length} \
            --input ~{input_bam} \
            --reference ~{reference_fa} \
            --tandem-repeats ~{tandems_bed} \
            --sample-id ~{sample_id} \
            --vcf ~{sample_id}_sniffles.vcf
        bgzip ~{sample_id}_sniffles.vcf
        tabix -f ~{sample_id}_sniffles.vcf.gz
    >>>

    output {
        File output_vcf_gz = sample_id + "_sniffles.vcf.gz"
        File output_tbi = sample_id + "_sniffles.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
