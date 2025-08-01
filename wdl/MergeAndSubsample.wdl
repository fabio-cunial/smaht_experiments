version 1.0


#
workflow MergeAndSubsample {
    input {
        Array[File] input_bams
        Array[File] input_bais
        Array[String] coverages
        String total_coverage
        String remote_output_dir
        String output_prefix
    }
    parameter_meta {
        coverages: "Comma-separated"
        total_coverage: "Precomputed estimate, to avoid running samtools coverage on the large merged file."
        output_prefix: "Of every file in `remote_output_dir`."
    }
    
    
    call Merge {
        input:
            input_bams = input_bams,
            input_bais = input_bais
    }
    scatter (i in range(length(coverages))) {
        call Subsample {
            input:
                merged_bam = Merge.merged_bam,
                merged_bai = Merge.merged_bai,
                coverage = coverages[i],
                total_coverage = total_coverage,
                remote_output_dir = remote_output_dir,
                output_prefix = output_prefix
        }
    }
    
    output {
    }
}


task Merge {
    input {
        Array[File] input_bams
        Array[File] input_bais
        
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
        
        INPUT_FILES=~{sep=',' input_bams}
        echo ${INPUT_FILES} | tr ',' '\n' > list.txt
        ${TIME_COMMAND} samtools merge --threads ${N_THREADS} -b list.txt -o merged.bam
        while read FILE; do
            rm -f ${FILE}
        done < list.txt
        ${TIME_COMMAND} samtools index --threads ${N_THREADS} merged.bam
    >>>
    
    output {
        File merged_bam = "merged.bam"
        File merged_bai = "merged.bam.bai"
    }
    runtime {
        docker: "fcunial/smaht_experiments"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


task Subsample {
    input {
        File merged_bam
        File merged_bai
        String coverage
        String total_coverage
        String remote_output_dir
        String output_prefix
        
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
        
        FRACTION=$(echo "scale=4; ~{coverage}/~{total_coverage}" | bc)
        ${TIME_COMMAND} samtools view --threads ${N_THREADS} --subsample ${FRACTION} --bam --output ~{output_prefix}_~{coverage}.bam ~{merged_bam}
        ${TIME_COMMAND} samtools index --threads ${N_THREADS} ~{output_prefix}_~{coverage}.bam
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m mv ~{output_prefix}_~{coverage}.'bam*' ~{remote_output_dir} && echo 0 || echo 1)
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
