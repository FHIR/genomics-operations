# Utilities
Utilities are used primarily to help load data into MongoDB and to support fast normalization.

## SPDI_Normalization
This code converts a chromosome-level variant, as derived from a VCF, into a contextual SPDI of the same build, using the algorithm described  [here](https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html). To run the code, you'll need to first download GRCh38 and GRCh37 Fasta files from  [NCBI Human Genome Resources page](https://www.ncbi.nlm.nih.gov/genome/guide/human/), and change the Python code to point to the downloaded files. The first time you run the code, the Fasta files get indexed, so it'll take longer.

## bed2json
Converts a BED file into a format suitable for loading into MongoDB. Chromosome numbering must include 'chr': 'chr1', 'chrX', 'chrY', 'chrM'. BED file must be sorted by chromosome, by position (bedtools sort default).

## vcf2json
Uses  [vcf2fhir](https://github.com/elimuinformatics/vcf2fhir)  logic to translate VCF records into a format suitable for loading into MongoDB.
