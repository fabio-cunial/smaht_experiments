version 1.0


#
workflow Gridss {
    input {
        String sample_id
        File input_bam
        File input_bai
        
        File reference_fa
        File reference_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
        
        Int n_cores = 16
        Int ram_gb = 32
    }
    parameter_meta {
    }

    call GridssImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reference_fa_amb = reference_fa_amb,
            reference_fa_ann = reference_fa_ann,
            reference_fa_bwt = reference_fa_bwt,
            reference_fa_pac = reference_fa_pac,
            reference_fa_sa = reference_fa_sa,
            n_cores = n_cores,
            ram_gb = ram_gb
    }

    output {
         File output_vcf_gz = GridssImpl.output_vcf_gz
         File output_tbi = GridssImpl.output_tbi
    }
}


#
task BwaIndex {
    input {
        File reference_fa
        File reference_fai
        
        Int n_cores = 4
        Int ram_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(reference_fa, "GB"))

    command <<<
        set -euxo pipefail
    
        TIME_COMMAND="/usr/bin/time --verbose"
        df -h
        
        mv ~{reference_fa} ./reference.fa
        mv ~{reference_fai} ./reference.fa.fai
        ${TIME_COMMAND} bwa index ./reference.fa
        ls -laht
    >>>

    output {
        File reference_fa_amb = "reference.fa.amb"
        File reference_fa_ann = "reference.fa.ann"
        File reference_fa_bwt = "reference.fa.bwt"
        File reference_fa_pac = "reference.fa.pac"
        File reference_fa_sa = "reference.fa.sa"
    }
    runtime {
        docker: "gridss/gridss:latest"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Remark: GRIDSS is run with default parameters. No crucial parameter seems to
# be exposed on the command line.
#
# Remark: we do not call `gridss_somatic_filter`, since it requires a panel of
# normals in a tumor/normal setting.
#
# Performance on a machine with 16 cores and 16 GB of RAM:
#
# COVERAGE  CPU%    TIME    RAM
# 120x
#
task GridssImpl {
    input {
        String sample_id
        File input_bam
        File input_bai
        
        File reference_fa
        File reference_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
        
        Int n_cores
        Int ram_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(input_bam, "GB")) + ceil(size(reference_fa, "GB")) + 200

    command <<<
        set -euxo pipefail
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        df -h
        
        ${TIME_COMMAND} gridss --threads ${N_THREADS} --jvmheap $(( ~{ram_gb} - 2 ))G --reference ~{reference_fa} --assembly assembly.bam --output out.vcf.gz ~{input_bam}
        ls -laht
        tree
        ${TIME_COMMAND} bcftools sort -m $(( ~{ram_gb} - 2 ))G --output-type z --output ~{sample_id}_gridss.vcf.gz out.vcf.gz
        ${TIME_COMMAND} tabix -f ~{sample_id}_gridss.vcf.gz
    >>>

    output {
        File output_vcf_gz = sample_id + "_gridss.vcf.gz"
        File output_tbi = sample_id + "_gridss.vcf.gz.tbi"
    }
    runtime {
        docker: "gridss/gridss:latest"
        cpu: n_cores
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
