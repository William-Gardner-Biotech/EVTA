# EVTA (Evolutionary Variant Timepoint Analysis)

[![Docker Image CI](https://github.com/William-Gardner-Biotech/EVTA/actions/workflows/docker-image.yml/badge.svg)](https://github.com/William-Gardner-Biotech/EVTA/actions/workflows/docker-image.yml)

EVTA is a bioinformatics pipeline for analyzing evolutionary changes in genomic variants across time intervals. It takes paired-end sequencing files (FASTQ) from two timepoints and uses the first timepoint as a reference to identify mutations that have occurred between the timespan.

## Overview

The pipeline consists of the following steps:

1. **Read QC and Merging**: The raw paired sequencing files (FASTQs) are quality-controlled and merged into a single file using bbduk and bbmerge.
2. **Map inoculum sample against SIVmac239M**: The inoculum sample (if provided, if not use the genbank reference) is mapped against the SIVmac239M reference genome to create an Inoculum Consensus FASTA.
3. **Map Timepoint 1 FASTQs to Inoculum reference**: The merged FASTQs from Timepoint 1 are mapped to the Inoculum Consensus FASTA to generate a Timepoint 1 Consensus FASTA.
4. **Map Timepoint 2 FASTQs to Timepoint 1 consensus**: The merged FASTQs from Timepoint 2 are mapped to the Timepoint 1 Consensus FASTA.
5. **Transfer SIVmac239M annotations to new consensus**: The annotations from the SIVmac239M reference are transferred to the Timepoint 2 BAM file using liftoff, creating a Revised GFF Annotation File.
6. **Call variants with the Revised Annotation File**: Variants are called using the Revised GFF Annotation File, generating a Variant Call Format (VCF) file.
Optional: Barcode region from SIVmac239M is removed using bedtools. This can be configured in the nextflow.config file.

## Usage

![Image 4-26-24 at 3 04 PM](https://github.com/William-Gardner-Biotech/EVTA/assets/99355149/f68bde09-d706-404f-abd5-e7132457ff23)


This repository contains sample data and sample output files to demonstrate the pipeline's functionality. To run the pipeline on your own data, follow these steps:

1. Clone the repository
2. Configure the Docker image
3. Place your paired-end sequencing files (FASTQs) in the `data/` directory
4. Mount volumes to Docker container
5. Configure your paths and settings in nextflow.config
6. Run the main.nf program

The pipeline will process your data, and the output files will be generated in the configured results directory.

## Docker image:

```
docker build -t evta .
```

OR

```
docker pull willgardnerbiotech/evta_1.0
```

## Contributing

Contributions are welcome! If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).
