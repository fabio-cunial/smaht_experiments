version 1.0


#
workflow Merge {
    input {
        Array[File] input_bams
        Array[File] input_bais

        String remote_output_dir
        String output_prefix
    }
    parameter_meta {
        output_prefix: "Of every file in `remote_output_dir`."
    }
    
    
    call Merge {
        input:
            input_bams = input_bams,
            input_bais = input_bais,
            remote_output_dir = remote_output_dir,
            output_prefix = output_prefix
    }
    
    output {
    }
}


# Performance on a set of 550x PacBio BAMs (~670GB total) on a VM with 8 cores
# and 32G of RAM:
#
# COMMAND               CPU%    TIME        RAM
# samtools merge        460%    2h46m       30M
# samtools index        200%    50m         200M
#
task Merge {
    input {
        Array[File] input_bams
        Array[File] input_bais
        
        String remote_output_dir
        String output_prefix
        
        Int n_cores = 6
        Int ram_gb = 4
        Int disk_size_gb = 1000
    }
    parameter_meta {
    }
    
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
        INPUT_FILES=~{sep=',' input_bams}
        echo ${INPUT_FILES} | tr ',' '\n' > list.txt
        ${TIME_COMMAND} samtools merge --threads ${N_THREADS} -b list.txt -o ~{output_prefix}.bam
        while read FILE; do
            rm -f ${FILE}
        done < list.txt
        ${TIME_COMMAND} samtools index --threads ${N_THREADS} ~{output_prefix}.bam
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m mv ~{output_prefix}.'bam*' ~{remote_output_dir} && echo 0 || echo 1)
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
