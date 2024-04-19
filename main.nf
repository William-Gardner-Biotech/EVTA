nextflow.enable.dsl=2

workflow {

    ch_fastqs_raw = Channel
        .fromFilePairs ( "${params.fastq_dir}/*{R1,R2}*.fastq.gz", size:2 )
        //.view()

    ch_samples = Channel
        .fromPath ( params.pair_file )
        .splitCsv ( header:true, strip:true )
        .map { row -> tuple(row.specimen, row.timepoint1, row.timepoint2) }
        .view()

    MERGE (
        ch_fastqs_raw
    )

    // Find the innoculum and build a consensus from it
    ASSEMBLE_INOC (
        MERGE.out
    )



    // merge all the paired reads together

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
    tuple val(sampleID), path("${sampleID}.fastq.gz")

    script:
    """
    bbmerge.sh in=${reads[0]} in2=${reads[1]} out=${sampleID}.fastq.gz
    """
}

// process to create a consensus from the positive control innoculum then return a consensus
process ASSEMBLE_INOC {

    tag "${sampleID}"

    publishDir "INOC", mode: 'symlink'

    input:
    tuple val(sampleID), path(merged_fastq)

    output:
    path "${params.innoc_name}.fa", emit: inoculum

    when:
    sampleID.contains(params.innoc_name)

    // We need to remove primers and trim adapters
    script:

    """
    bbduk.sh -Xmx8g in=${merged_fastq} out=QC_${merged_fastq} minlen=100 hdist=2 ftm=5 maq=10
    bbmap.sh -Xmx8g ref=${params.ref_fasta} in=QC_${merged_fastq} out=${sampleID}.bam maxindel=100 minid=0.9
    samtools sort -o sorted_${sampleID}.bam -l 1 ${sampleID}.bam
    samtools mpileup sorted_${sampleID}.bam | ivar consensus -p ${params.innoc_name}
    """

}

