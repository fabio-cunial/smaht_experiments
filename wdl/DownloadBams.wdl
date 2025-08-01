version 1.0


#
workflow DownloadBams {
    input {
        File access_key_csv
        File manifest_tsv
        String remote_output_dir
        
        Int n_cores
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
        manifest_tsv: "From the data portal"
        remote_dir: "Output directory in a remote bucket"
    }
    
    call DownloadBamsImpl {
        input:
            access_key_csv = access_key_csv,
            manifest_tsv = manifest_tsv,
            remote_output_dir = remote_output_dir,
            n_cores = n_cores,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
        File coverage_txt = DownloadBamsImpl.coverage_txt
    }
}


task DownloadBamsImpl {
    input {
        File access_key_csv
        File manifest_tsv
        String remote_output_dir
        
        Int n_cores
        Int ram_gb
        Int disk_size_gb
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
        
        # Retrieving access key
        ACCESS_KEY_ID=$(head -n 1 ~{access_key_csv} | cut -d , -f 1)
        ACCESS_KEY_SECRET=$(head -n 1 ~{access_key_csv} | cut -d , -f 2)
        
        # Downloading from the manifest
        cut -f 1,3 ~{manifest_tsv} | tail -n +4 | tr '\t' ',' > list.txt
        while read LINE; do
            ADDRESS=$(echo ${LINE} | cut -d , -f 1)
            FILENAME=$(echo ${LINE} | cut -d , -f 2)
            curl -L --user ${ACCESS_KEY_ID}:${ACCESS_KEY_SECRET} ${ADDRESS} --output ${FILENAME} 
        done < list.txt
        ${TIME_COMMAND} samtools coverage *.bam > coverage.txt
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp '*.bam*' ~{remote_output_dir} && echo 0 || echo 1)
            if [[ ${TEST} -eq 1 ]]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
        File coverage_txt = "coverage.txt"
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
