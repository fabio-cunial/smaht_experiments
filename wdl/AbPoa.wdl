version 1.0


#
workflow AbPoa {
    input {
        File reads_fa
        Int abpoa_align_mode = 1
        Int abpoa_result = 3
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            reads_fa = reads_fa,
            abpoa_align_mode = abpoa_align_mode,
            abpoa_result = abpoa_result
    }
    
    output {
        File out_txt = Impl.out_txt
    }
}


task Impl {
    input {
        File reads_fa
        Int abpoa_result
        Int abpoa_align_mode
        
        Int n_cores = 2
        Int ram_gb = 128
    }
    parameter_meta {
        abpoa_align_mode: "The `--aln-mode` flag of abPOA: 0=global, 1=local, 2=extension."
        abpoa_result: "The `--result` flag of abPOA: 1=MSA 3=GFA."
    }
    
    Int disk_size_gb = 10*ceil(size(reads_fa, "GB"))
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
        
        ${TIME_COMMAND} ~{docker_dir}/abPOA/bin/abpoa -m ~{abpoa_align_mode} --amb-strand --sort-by-len --result ~{abpoa_result} --verbose 2 ~{reads_fa} > out.txt
        ls -laht
    >>>
    
    output {
        File out_txt = "out.txt"
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
