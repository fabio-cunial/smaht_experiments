version 1.0


#
workflow Severus {
    input {
        File aligned_bam
        File aligned_bai
        
        File pon_tsv_gz
        File vntr_bed
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            pon_tsv_gz = pon_tsv_gz,
            vntr_bed = vntr_bed
    }
    
    output {
        File out_tar_gz = Impl.out_tar_gz
    }
}


# 
#
# Performance on a machine with 32 cores and 64GB of RAM:
#
# COVERAGE  CPU%    TIME    RAM
# 230x      760%    1h40m   25G
#
task Impl {
    input {
        File aligned_bam
        File aligned_bai
        
        File pon_tsv_gz
        File vntr_bed
        
        Int n_cores = 8
        Int ram_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(aligned_bam, "GB") + size(pon_tsv_gz, "GB")) 
    String docker_dir = "/smaht_experiments"
    
    command <<<
        set -euxo pipefail
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        df -h
        
        source activate severus_env
        ${TIME_COMMAND} severus --threads ${N_THREADS} --target-bam ~{aligned_bam} --PON ~{pon_tsv_gz} --vntr-bed ~{vntr_bed} --out-dir ./severus 
        ls -laht
        ${TIME_COMMAND} tar -czf out.tar.gz ./severus
    >>>
    
    output {
        File out_tar_gz = "out.tar.gz"
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
