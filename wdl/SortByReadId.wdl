version 1.0


#
workflow SortByReadId {
    input {
        File input_bam
        File input_bai
        
        String remote_output_dir
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            input_bam = input_bam,
            input_bai = input_bai,
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}


# Performance on a 230x PacBio BAMs on a VM with 6 cores and 16G of RAM:
#
# COMMAND               CPU%    TIME        RAM
# 
#
task Impl {
    input {
        File input_bam
        File input_bai
        
        String remote_output_dir
        
        Int n_cores = 6
        Int ram_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_bam, "GB"))
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
        
        # Merging
        FILENAME=$(basename ~{input_bam})
        FILENAME=${FILENAME%.*}
        ${TIME_COMMAND} samtools collate -@ ${N_THREADS} -o ${FILENAME}_sorted.bam ~{input_bam}
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ${FILENAME}_sorted.bam
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m mv ${FILENAME}_sorted.'bam*' ~{remote_output_dir} && echo 0 || echo 1)
            if [[ ${TEST} -eq 1 ]]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
