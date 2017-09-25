##############################################################################
# Workflow Definition
##############################################################################

workflow EpitopePrediction {
    File inputs_tsv
    Array[Array[File]] read_pairs_array = read_tsv(inputs_tsv)
    String output_basename
    
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index

    scatter (read_pairs in read_pairs_array) {
        call GetRgids {
            input:
                first = read_pairs[0],
                output_basename = output_basename
        }
        call BwaMemAndSortFixTags {
            input:
                first = read_pairs[0],
                second = read_pairs[1],
                rgids = GetRgids.rgids,
                output_bam_basename = output_basename,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                ref_bwt = ref_bwt,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_pac = ref_pac,
                ref_sa = ref_sa
        }
    }
    call MergeBamFiles {
        input:
            input_bams = BwaMemAndSortFixTags.output_bam,
            output_bam_basename = output_basename
    }
    call BaseRecalibrator {
        input:
            input_bam = MergeBamFiles.output_bam,
            input_bam_index = MergeBamFiles.output_bam_index,
            recalibration_report_filename = output_basename
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            known_indels = known_indels,
            known_indels_index = known_indels_index,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }
    call ApplyBQSR {
        input:
            input_bam = MergeBamFiles.output_bam,
            input_bam_index = MergeBamFiles.output_bam_index,
            output_bam_basename = output_basename,
            recalibration_report = BaseRecalibrator.recalibration_report,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
    }
    call HaplotypeCaller {
        input:
            input_bam = ApplyBQSR.recalibrated_bam,
            input_bam_index = ApplyBQSR.recalibrated_bam_index,
            gvcf_name = output_basename,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }
}

##############################################################################
# Task Definitions
##############################################################################

task GetRgids {
    File first
    String output_basename

    command {
        python ~/bin/get_RGID_from_FASTQ.py first output_basename
    }
    output {
        String rgids = read_string(stdout())
    }
}

task BwaMemAndSortFixTags {
    File first
    File second
    String rgids
    String output_bam_basename

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    command {
        ~/bin/bwa mem -p -v 3 -t 3 -R ${rgids} \
            ${ref_fasta} \
            ${first} ${second} \
        2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
        > ${output_bam_basename}.sam;
        ~/bin/java -Xmx2500m -jar ~/bin/picard.jar \
            SamFormatConverter \
            I=${output_bam_basename}.sam \
            O=/dev/stdout \
        ~/bin/java -Xmx8000m -jar ~/bin/picard.jar \
            SortSam \
            INPUT=/dev/stdin \
            OUTPUT=/dev/stdout \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false | \
        java -Xmx1000m -jar /usr/bin_dir/picard.jar \
            SetNmAndUqTags \
            INPUT=/dev/stdin \
            OUTPUT=${output_bam_basename}.sorted.bam \
            CREATE_INDEX=true \
            CREATE_MD5_FILE=false \
            REFERENCE_SEQUENCE=${ref_fasta};
    }
    output {
        File output_bam = "${output_bam_basename}.sorted.bam"
        File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
    }
}

task MergeBamFiles {
    Array[File] input_bams
    String output_bam_basename

    command {
        ~/bin/java -Xmx8000m -jar ~/bin/picard.jar \
            MergeSamFiles \
            I=${sep " I=" input_bams} \
            O=${output_bam_basename}.bam \
            SORT_ORDER=coordinate \
            ASSUME_SORTED=true;
    }
    output {
        File output_bam = "${output_bam_basename}.bam"
        File output_bam_index = "${output_bam_basename}.bai"
    }
}

task BaseRecalibrator {
    File input_bam
    File input_bam_index
    String recal_report_filename
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int disk_size
    Int preemptible_tries

    command {
        java -Xmx4000m -jar /usr/bin_dir/GATK.jar \
            -T BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -o ${recal_report_filename}.recal_data.table \
            -knownSites ${dbsnp} \
            -knownSites ${known_indels}
    }
    output {
        File recalibration_report = "${recal_report_filename}.recal_data.table"
    }
}

task ApplyBQSR {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int disk_size
    Int preemptible_tries

    command {
        java -Xmx2000m -jar /usr/bin_dir/GATK.jar \
            -T PrintReads \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -o ${output_bam_basename}.realn.sorted.bqsr.bam \
            -BQSR ${recalibration_report}
    }
    output {
        File recalibrated_bam = "${output_bam_basename}.sorted.bqsr.bam"
        File recalibrated_bam_index = "${output_bam_basename}.sorted.bqsr.bai"
    }
}

task HaplotypeCaller {
    File input_bam
    File input_bam_index
    String vcf_name
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    command {
        java -Xmx8000m -jar /usr/bin_dir/GATK.jar \
            -T HaplotypeCaller \
            -R ${ref_fasta} \
            -o ${vcf_name}.vcf \
            -I ${input_bam}
    }
    output {
        File output_gvcf = "${vcf_name}.vcf"
        #File output_gvcf_index = "${vcf_name}.vcf.tbi"
    }
}
