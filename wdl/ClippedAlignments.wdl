version 1.0


#
workflow ClippedAlignments {
    input {
        File input_bam
        
        File script_java
        String remote_output_dir
    }
    parameter_meta {
        input_bam: "Sorted by read ID"
    }
    
    
    call Impl {
        input:
            input_bam = input_bam,
            script_java = script_java,
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}


#
task Impl {
    input {
        File input_bam
        
        File script_java
        String remote_output_dir
        
        Int n_cores = 2
        Int ram_gb = 16
    }
    parameter_meta {
        input_bam: "Sorted by read ID"
        remote_output_dir: "Must end with a '/', otherwise it is interpreted as a file if it does not exist."
    }
    
    Int disk_size_gb = 5*( ceil(size(input_bam,"GB")) )
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
        
        # Printing clipped reads
        FILE_NAME=$(basename ~{input_bam})
        mv ~{script_java} ./ClippedAlignments.java
        javac ClippedAlignments.java
        date
        samtools view ~{input_bam} | java ClippedAlignments > ${FILE_NAME}_clipped_reads.csv
        date
        ${TIME_COMMAND} gzip ${FILE_NAME}_clipped_reads.csv
        ls -laht
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m mv '*_clipped_reads.csv.gz' ~{remote_output_dir} && echo 0 || echo 1)
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
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
