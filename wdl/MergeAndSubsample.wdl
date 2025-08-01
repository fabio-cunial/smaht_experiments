version 1.0


#
workflow MergeAndSubsample {
    input {
        Array[File] input_bams
        Array[File] input_bais
        String coverages
        String total_coverage
        String remote_output_dir
        String output_prefix
        
        Int n_cores
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
        coverages: "Comma-separated"
        total_coverage: "Precomputed estimate, to avoid running samtools coverage on the large merged file."
        output_prefix: "Of every file in `remote_output_dir`."
        n_cores: "Exactly this number of threads will be used."
    }
    
    call MergeAndSubsampleImpl {
        input:
            input_bams = input_bams,
            input_bais = input_bais,
            coverages = coverages,
            total_coverage = total_coverage,
            n_cores = n_cores,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


task MergeAndSubsampleImpl {
    input {
        Array[File] input_bams
        Array[File] input_bais
        String coverages
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
        df -h
        
        function subsample() {
            local MERGED_BAM=$1
            local COVERAGES_FILE=$2
            
            while read COVERAGE; do
                FRACTION=$(echo "scale=4; ${COVERAGE}/~{total_coverage}" | bc)
                ${TIME_COMMAND} samtools view --threads 1 --subsample ${FRACTION} --bam --output ~{output_prefix}_${COVERAGE}.bam ${MERGED_BAM}
                ${TIME_COMMAND} samtools index --threads 1 ~{output_prefix}_${COVERAGE}.bam
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m mv ~{output_prefix}_${COVERAGE}.'bam*' ~{remote_output_dir} && echo 0 || echo 1)
                    if [[ ${TEST} -eq 1 ]]; then
                        echo "Error uploading files. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            done < ${COVERAGES_FILE}
        }
        
        # Merging
        echo "${sep='\n' input_bams}" > list.txt
        ${TIME_COMMAND} samtools merge --threads ${N_THREADS} -b list.txt -o merged.bam
        while read FILE; do
            rm -f ${FILE}
        done < list.txt
        ${TIME_COMMAND} samtools index --threads ${N_THREADS} merged.bam
        
        # Subsampling
        echo ~{coverages} | tr ',' '\n' | sort --random-sort > list.txt
        N_COVERAGES=$(wc -l < list.txt)
        N_COVERAGES_PER_THREAD=$( ${N_COVERAGES} / ~{n_cores} )
        split -d -l ${N_COVERAGES_PER_THREAD} list.txt chunk_
        for FILE in $(ls chunk_); do
            subsample merged.bam ${FILE} &
        done
        wait
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
