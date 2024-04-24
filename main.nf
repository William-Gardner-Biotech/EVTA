nextflow.enable.dsl=2

log.info """
                                                                      
8 8888888888 `8.`888b           ,8' 8888888 8888888888   .8.          
8 8888        `8.`888b         ,8'        8 8888        .888.         
8 8888         `8.`888b       ,8'         8 8888       :88888.        
8 8888          `8.`888b     ,8'          8 8888      . `88888.       
8 888888888888   `8.`888b   ,8'           8 8888     .8. `88888.      
8 8888            `8.`888b ,8'            8 8888    .8`8. `88888.     
8 8888             `8.`888b8'             8 8888   .8' `8. `88888.    
8 8888              `8.`888'              8 8888  .8'   `8. `88888.   
8 8888               `8.`8'               8 8888 .888888888. `88888.  
8 888888888888        `8.`                8 8888.8'       `8. `88888. 


""".stripIndent()

workflow {

    ch_fastqs_raw = Channel
        .fromFilePairs ( "${params.fastq_dir}/*{R1,R2}*.fastq.gz", size:2 )
        //.view()

    ch_samples = Channel
        .fromPath ( params.pair_file )
        .splitCsv ( header:true, strip:true )
        .map { row -> tuple(row.specimen, row.timepoint1, row.timepoint2) }
        //.view()

    MERGE (
        ch_fastqs_raw
    )

    // Find the innoculum and build a consensus from it
    ASSEMBLE_INOC (
        MERGE.out
    )

    // Working channel that combines merged samples by the specimen ID key
    // Defined by the everything in the sample name before a dash char (like .split('-')[0])
    grouped_ch = MERGE.out
        .groupTuple(
            size: 2)
            .map { x -> tuple(x[0], x[1][0], x[1][1]) } //unpack list
            //.view()

    GENERATE_T1_CONSENSUS (
        grouped_ch,
        ASSEMBLE_INOC.out
    )

    LIFT_GFF (
        GENERATE_T1_CONSENSUS.out
    )

    if (params.snp_eff_db_exists == false) {
            MAKE_SNPEFF_DB (
                LIFT_GFF.out
                )
    }

    MAP_SECOND_TIMEPOINT (
        LIFT_GFF.out
    )

    CALL_VARIANTS (
        MAP_SECOND_TIMEPOINT.out
    )



    // USE Annotation Transfer Tool to transfer annotation to new consensus genomes

    // Make snpEff database to reference
    // Use the snpEff build -noCheckCds -noCheckProtein

    // Call variants using SNPeff
    
    // Count synonymous vs non and return the data

    // QC the reads with decent strictness before variant calling


}

// MERGE the paired reads together and check to find the inoculum to emit it elsewhere
// Also remove Nextera adapters
process MERGE {

    tag "${sampleID}"

    publishDir "Merged", mode: 'copy'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sample), path("${sampleID}.fastq.gz")

    script:

    sample = sampleID.tokenize('-')[0]
    """
    bbmerge.sh in=${reads[0]} in2=${reads[1]} out=${sampleID}.fastq.gz
    """
}

// process to create a consensus from the positive control innoculum then return a consensus
process ASSEMBLE_INOC {

    tag "${sampleID}"

    publishDir "INOC", mode: 'symlink'

    input:
    tuple val(specimen), path(merged_fastq)

    output:
    path "${params.inoc_prefix_before_dash}.fa"

    when:
    specimen.contains(params.inoc_prefix_before_dash)

    // We need to remove primers and trim adapters
    script:

    """
    bbduk.sh -Xmx8g in=${merged_fastq} out=QC_${merged_fastq} minlen=100 hdist=2 ftm=5 maq=10
    bbmap.sh -Xmx8g ref=${params.ref_fasta} in=QC_${merged_fastq} out=${specimen}.bam maxindel=100 minid=0.9
    samtools sort -o sorted_${specimen}.bam -l 1 ${specimen}.bam
    samtools mpileup sorted_${specimen}.bam | ivar consensus -p ${params.inoc_prefix_before_dash}
    """

}

// This process will take the timepoint pair channel and then match the corresponding merged reads with it
process CROSS_CHANNELS {
    input:
    tuple val(specimen), val(T1), val(T2)
    tuple val(sampleID), path(merged_fastq)

    output:
    tuple val(specimen), path(matched_read1), path(matched_read2)

    script:
    read1_name = T1
    read2_name = T2

    matched_read1 = merged_fastq.findAll { it.getName()
    .startsWith(read1_name)
    }
    matched_read2 = merged_fastq.findAll { it.getName()
    .startsWith(read2_name)
    }

    if (matched_read1.size() == 1 && matched_read2.size() == 1) {
        """
        echo "${specimen} $matched_read1 $matched_read2"
        """
    } else {
        println "Skipping ${specimen} as one or both files could not be matched uniquely."
    }

}

process GENERATE_T1_CONSENSUS {

    tag "${specimen}"

    publishDir "T1_Consensus", mode: 'Copy'

    // T = Timepoint1
    input:
    tuple val(specimen), path(T1), path(T2)
    each path(inocolum)

    output:
    tuple val(specimen), path(timepoint2), path("${specimen}.fa")

    script:

    // if the sample number for T1 is larger than T2 then we need to swap the vars for timewise comparison
    if (T1.getSimpleName().tokenize('-')[1] > T2.getSimpleName().tokenize('-')[1]) {
        timepoint1 = T2
        timepoint2 = T1
    }
    else {
        timepoint1 = T1
        timepoint2 = T2
    }

    """
    echo "${timepoint1} | ${timepoint2}"
    bbduk.sh -Xmx8g in=${timepoint1} out=QC_${timepoint1} minlen=100 hdist=2 ftm=5 maq=10
    bbmap.sh -Xmx8g ref=${params.ref_fasta} in=QC_${timepoint1} out=${specimen}.bam maxindel=100 minid=0.9
    samtools sort -o sorted_${specimen}.bam -l 1 ${specimen}.bam
    samtools mpileup sorted_${specimen}.bam | ivar consensus -p ${specimen}
    """
}

// process that will take ref GFF and ref fasta then compare new consensus to shift the annotations
process LIFT_GFF {

    publishDir "GFF_lift", mode: 'Copy'

    tag "${specimen}"

    input:
    tuple val(specimen), path(timepoint_2), path(t1_Consensus)

    output:
    tuple val(specimen), path(timepoint_2), path(t1_Consensus), path("${specimen}.gff3")

    script:

    """
    liftoff -g ${params.ref_GFF} -o ${specimen}.gff3 ${t1_Consensus} ${params.ref_fasta}
    """
}

// process that needs to be run once for the sake of making a snpEff database to call upon
// If you run it multiple times it will have collision errors
// TODO make this smarter
process MAKE_SNPEFF_DB {

    tag "${specimen}"

    input:
    tuple val(specimen), path(timepoint_2), path(t1_Consensus), path(consensus_gff)

    output:
    val specimen

    script:

    """

    # Make the database now
    echo "\n# ${specimen} genome\n${specimen}.genome: ${specimen}" >> ${params.snpEff_folder}/snpEff.config
    mkdir ${params.snpEff_folder}/data/${specimen}
    cp ${t1_Consensus} ${params.snpEff_folder}/data/${specimen}/sequences.fa
    cp ${specimen}.gff3 ${params.snpEff_folder}/data/${specimen}/genes.gff

    # Build the database with snpEff
    java -jar ${params.snpEff_folder}/snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v ${specimen}
    """
}

// Also mapping the amp3 primers we care about to use bcftools subtract
process MAP_SECOND_TIMEPOINT {

    tag "${specimen}"

    publishDir "SecondBam", mode: 'Copy'

    input:
    tuple val(specimen), path(timepoint_2), path(t1_Consensus), path(consensus_gff)

    output:
    tuple val(specimen), path("${specimen}_vcf_callable.bam"), path(t1_Consensus), path(consensus_gff)

    // Mapping portion that also removes the amplicon regions not represented by all samples (specified in amp_primers.fa)
    // Then it maps and removes the barcode region from the bam file and converts all the way back to a final bam file.
    script:
    """
    # Map the amp primers so we can subtract down full bam
    bbmap.sh ref=${t1_Consensus} in=${params.amp_primers} out=mapped_primers.bam
    bedtools bamtobed -i mapped_primers.bam > mapped_primers.bed

    awk 'BEGIN {OFS="\t"} NR==1 {chr=\$1; start=\$2; name1=\$4; score1=\$5; strand=\$6} END {print chr "\t" start "\t" \$3 "\t" name1 "\t" score1 "\t" strand}' mapped_primers.bed > primer_amplicon.bed

    # Map the second timepoint reads against time1 consensus
    bbduk.sh -Xmx8g in=${timepoint_2} out=QC_${timepoint_2} minlen=100 hdist=2 ftm=5 maq=10
    bbmap.sh -Xmx8g ref=${t1_Consensus} in=QC_${timepoint_2} out=${specimen}_full.bam maxindel=100 minid=0.9

    # Subtract the reads outside the primer region
    # bedtools bamtobed -i ${specimen}_full.bam > ${specimen}_full.bed
    ###################################################################################################
    # bedtools intersect -a ${specimen}_full.bam -b primer_amplicon.bed > ${specimen}_shrunk.bam
    # ################### ^^ Line is important for filtering down to single amplicon ^^ ###############

    # Remove the barcode region
    awk '\$3 == "gene" && (\$9 ~ /vpr/ || \$9 ~ /vpx/)' ${consensus_gff} | awk '{print \$1, \$4, \$5, \$7}' > vpr_vpx_coords.txt

    # awk will then search for vpr-vpx gene region and print it into a bed file.
    # It also adds a 15 bp flank to each side to remove some mapping errors near the region
    awk 'NR==1 {vpr_start=\$2; vpr_end=\$3; vpr_strand=\$4; next} 
        NR==2 {vpx_start=\$2; vpx_end=\$3; vpx_strand=\$4; 
                if (vpx_strand == "+") {
                    start = vpr_end + 1 - 15;
                    end = vpx_start - 1 + 15;
                    print \$1 "\t" start "\t" end "\t" vpx_strand
                } else {
                    start = vpx_end + 1 - 15;
                    end = vpr_start - 1 + 15;
                    print \$1 "\t" start "\t" end "\t" vpr_strand
                }
        }' vpr_vpx_coords.txt > between_vpr_vpx.bed

    bedtools subtract -a ${specimen}_full.bam -b between_vpr_vpx.bed -A > ${specimen}_vcf_callable.bam 
    """
}

process CALL_VARIANTS {

    tag "${specimen}"

    input:
    tuple val(specimen), path(callable_bam), path(t1_Consensus), path(consensus_gff)

    output:
    tuple val(specimen), path("${specimen}.tsv"), path(t1_Consensus), path(consensus_gff)

    publishDir "Final_tsv_SNPs", mode: 'Copy'

    script:
    """
    samtools sort -o sorted_${callable_bam} -l 1 ${callable_bam}
    samtools mpileup sorted_${callable_bam} | ivar variants -r ${t1_Consensus} -p ${specimen} -g ${consensus_gff} -m 10
    """
}

/*
We need to subtract the barcode region from the bam file, and 15 bp upstream and down from that. We can use the gff file to find it as it lies
betweem vpx and vpr genes. 
cy0661 only got amplicon 3 so we need to filter all vcf's down to that region.
*/