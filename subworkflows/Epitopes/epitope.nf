nextflow.enable.dsl=2

log.info """
                                                                    
                                                                                     
8 8888888888 `8.`888b           ,8' 8888888 8888888888   .8.                8 8888888888   
8 8888        `8.`888b         ,8'        8 8888        .888.               8 8888         
8 8888         `8.`888b       ,8'         8 8888       :88888.              8 8888         
8 8888          `8.`888b     ,8'          8 8888      . `88888.             8 8888         
8 888888888888   `8.`888b   ,8'           8 8888     .8. `88888.            8 888888888888 
8 8888            `8.`888b ,8'            8 8888    .8`8. `88888.           8 8888         
8 8888             `8.`888b8'             8 8888   .8' `8. `88888.          8 8888         
8 8888              `8.`888'              8 8888  .8'   `8. `88888.         8 8888         
8 8888               `8.`8'               8 8888 .888888888. `88888.   d8b  8 8888         
8 888888888888        `8.`                8 8888.8'       `8. `88888.  Y8P  8 888888888888 


Epitopes Only

${params.ref_fasta}

""".stripIndent()

// Before your workflow, add the CSV reading function:
def readThresholds(csv_file) {
    def thresholds = [:]
    new File(csv_file).withReader { reader ->
        def header = reader.readLine()  // Skip header
        reader.eachLine { line ->
            def (specimen, threshold) = line.split(',')
            thresholds[specimen] = threshold
        }
    }
    return thresholds
}

workflow {

    ch_fastqs_raw = Channel
        .fromFilePairs ( "${params.fastq_dir}/*{R1,R2}*.fastq.gz", size:2 )
        //.view()

    MERGE (
        ch_fastqs_raw
    )

    MAP_READS_ON_SIV (
        MERGE.out,
        params.ref_fasta,
        params.ref_GFF
    )

    EXTRACT_EPITOPES (
        params.epitope_gff
    )

    //Read CSV
    thresholds = readThresholds(params.vcf_csv)

    // Create a new channel with BAM files and their corresponding thresholds
    ch_bam_with_threshold = MAP_READS_ON_SIV.out.map { specimen, bam ->
        // thresholds is a dictionary which we are now indexing
        def threshold = thresholds[specimen] ?: params.freq_th  // Default to frequency threshold if not found
        tuple(specimen, bam, threshold)
    }

    SPECIALIZED_VARIANT_CALL (
        ch_bam_with_threshold,
        EXTRACT_EPITOPES.out
    )

}

// MERGE the paired reads together
process MERGE {

    tag "${sampleID}"

    publishDir "Merged", mode: 'copy'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sample), path("${sample}.fastq.gz")

    script:

    sample = sampleID.tokenize('-')[0]
    """
    bbmerge.sh in=${reads[0]} in2=${reads[1]} out=${sample}.fastq.gz
    """
}

// Also mapping the amp3 primers we care about to use bcftools subtract
process MAP_READS_ON_SIV {

    tag "${specimen}"

    publishDir "Sorted_n_Mapped_BAMs", mode: 'Copy'

    input:
    tuple val(specimen), path(merged_reads)
    each path(ref_seq)
    each path(ref_gff)

    output:
    tuple val(specimen), path("sorted_${specimen}.bam")

    // Mapping portion that also removes the amplicon regions not represented by all samples (specified in amp_primers.fa)
    // Then it maps and removes the barcode region from the bam file and converts all the way back to a final bam file.
    script:
    """

    # QC then map the reads onto our reference
    bbduk.sh -Xmx8g in=${merged_reads} out=QC_${merged_reads} minlen=100 hdist=2 ftm=5 maq=10
    bbmap.sh -Xmx8g ref=${ref_seq} in=QC_${merged_reads} out=${specimen}_mapped.bam maxindel=100 minid=0.9

    # Remove the barcode region
    awk '\$3 == "gene" && (\$9 ~ /vpr/ || \$9 ~ /vpx/)' ${ref_gff} | awk '{print \$1, \$4, \$5, \$7}' > vpr_vpx_coords.txt

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

    bedtools subtract -a ${specimen}_mapped.bam -b between_vpr_vpx.bed -A > ${specimen}_mapped_no_bc.bam 

    samtools sort -o sorted_${specimen}.bam -l 1 ${specimen}_mapped_no_bc.bam
    #samtools mpileup sorted_${specimen}.bam | ivar variants -r ${ref_seq} -p ${specimen}_t1 -g ${ref_gff} -m 100
    """
}

// TODO make this automated
process EXTRACT_EPITOPES {

    publishDir "Epitopes_only_gff", mode: 'symlink'

    input:
    path gff_file

    output:
    path "epitopes.bed"

    script:
    """
    awk -F'\t' '\$3 == "misc_feature" {
        split(\$9, attrs, ";");
        for (i in attrs) {
            if (attrs[i] ~ /^Name=/) {
                name = substr(attrs[i], 6);
                break;
            }
        }
        print \$1 "\t" \$4-1 "\t" \$5 "\t" name "\t" "+"
    }' ${gff_file} > epitopes.bed
    """
}

process SPECIALIZED_VARIANT_CALL {

    tag "Variant calling: ${specimen} using a -t=${threshold}"

    publishDir "Variant_calls", mode: 'Copy'

    input:
    tuple val(specimen), path(sorted_bam), val(threshold)
    each path(epitope_bed)

    output:
    tuple val(specimen), path("${specimen}.tsv")

    script:
    """
    # ivar -m is min depth, -t is variant threshold, -r ref seq, -p is prefix for output file
    samtools mpileup -f ${params.ref_fasta} ${sorted_bam} | ivar variants -r ${params.ref_fasta} -p ${specimen} -g ${params.ref_GFF} -m 100 -t ${threshold}

    """
}