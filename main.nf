nextflow.enable.dsl=2

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
    tuple val(specimen), path(T2), path("${specimen}.fa")

    // if the sample number for T1 is larger than T2 then we need to swap the vars for timewise comparison
    if (T1.getSimpleName().tokenize('-')[1] > T2.getSimpleName().tokenize('-')[1]) {
        def temp = T1
        T1 = T2
        T2 = temp
    }

    """
    echo ${T1} | ${T2}
    bbduk.sh -Xmx8g in=${T1} out=QC_${T1} minlen=100 hdist=2 ftm=5 maq=10
    bbmap.sh -Xmx8g ref=${params.ref_fasta} in=QC_${T1} out=${specimen}.bam maxindel=100 minid=0.9
    samtools sort -o sorted_${specimen}.bam -l 1 ${specimen}.bam
    samtools mpileup sorted_${specimen}.bam | ivar consensus -p ${specimen}
    """
}

process BUILD_snpEff_DB {

    tag "${specimen}"

    input:
    tuple val(specimen), path(T2), path(t1_Consensus)

    output:
    tuple val(specimen), path(snpEff_db)

    script:

    """
    liftoff -g ${params.ref_GFF} -o ${specimen}.gff3 ${t1_Consensus} ${params.ref_fasta}

    # Make the database now
    echo "# ${specimen} genome\n${specimen}.genome: ${specimen}" >> ${params.snpEff_folder}/snpEff.config
    mkdir ${params.snpEff_folder}/data/${specimen}
    cp ${t1_Consensus} ${params.snpEff_folder}/data/${specimen}/sequences.fa
    cp ${specimen}.gff3 ${params.snpEff_folder}/data/${specimen}/genes.gff

    # Build the database with snpEff
    java -jar ${params.snpEff_folder}/snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v ${specimen}
    """
}
