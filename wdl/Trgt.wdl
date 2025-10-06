version 1.0


#
workflow Trgt {
    input {
        File aligned_bam
        File aligned_bai
        File reference_fa
        File reference_fai
        File catalog_bed
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            catalog_bed = catalog_bed
    }
    
    output {
        Array[File] out = Impl.out
    }
}


# Performance on a machine with 16 cores and 32GB of RAM:
#
# COVERAGE  CPU%    TIME    RAM
# 230x      
#
task Impl {
    input {
        File aligned_bam
        File aligned_bai
        File reference_fa
        File reference_fai
        File catalog_bed
        
        Int n_cores = 16
        Int ram_gb = 32
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
        
        
        ${TIME_COMMAND} ~{docker_dir}/trgt genotype --threads ${N_THREADS} --genome ~{reference_fa} --reads ~{aligned_bam} --repeats ~{catalog_bed} --output-prefix out 
        ls -laht
    >>>
    
    output {
        Array[File] out = glob("out*")
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
