version 1.0


#
workflow Savana {
    input {
        File aligned_bam
        File aligned_bai
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    
    output {
        File out_tar_gz = Impl.out_tar_gz
    }
}


# Remark: savana already claims to use the maximum possible number of theads.
#
# Performance on a machine with 32 cores and 64GB of RAM:
#
# COVERAGE  CPU%    TIME    RAM
# 230x      2800%   >=6h    20G
#
task Impl {
    input {
        File aligned_bam
        File aligned_bai
        File reference_fa
        File reference_fai
        
        Int n_cores = 32
        Int ram_gb = 64
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(aligned_bam, "GB") + size(reference_fa, "GB")) 
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
        
        
        ${TIME_COMMAND} savana to --pb --tumour ~{aligned_bam} --ref ~{reference_fa} --g1000_vcf 1000g_hg38 --contigs ~{docker_dir}/savana/example/contigs.chr.hg38.txt --outdir ./savana --sample out
        ls -laht
        ${TIME_COMMAND} tar -czf out.tar.gz ./savana
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
