# Utilities
Utilities are used primarily to help load data into MongoDB and to support fast normalization.

## SPDI_Normalization
This code converts a chromosome-level variant, as derived from a VCF, into a contextual SPDI of the same build, using the algorithm described  [here](https://web.archive.org/web/20240814131742/https://vrs.ga4gh.org/en/stable/impl-guide/normalization.html). To run the code, you'll need to first download GRCh38 and GRCh37 Fasta files (from [this page](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/), navigate to latest GRCh37 and GRCh38 assemblies (e.g. [GRCh37](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/); [GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/)). For each assembly, download the FASTA file, which is named "..._genomic.fna.gz") and change the Python code to point to the downloaded files. The first time you run the code, the Fasta files get indexed, so it'll take longer.

## bed2json
Converts a BED file into a format suitable for loading into MongoDB. Chromosome numbering must include 'chr': 'chr1', 'chrX', 'chrY', 'chrM'. BED file must be sorted by chromosome, by position (bedtools sort default).

## run_vcf2json
Batch process that calls vcf2json for a set of VCF files, yielding three output files ('variantsData.json', 'phaseData.json', 'molecularConsequences.json') for loading into respective MongoDB collections. Does not update Patients or Tests collections. VCFs to be processed are listed in vcfData.csv, which must include columns _vcf_filename_, _ref_build_ (populated with 'GRCh37' or 'GRCh38'), _patient_id_, _test_date_ (yyyy-mm-dd), _test_id_, _specimen_id_, _genomic_source_class_ (populate with 'germline', 'somatic', or 'mixed'), _ratio_ad_dp_ (used for mitochondrial DNA processing, generally set it to 0.99), _sample_position_ (zero-based, useful for multi-sample VCFs). vcf2json translation logic is based on [vcf2fhir](https://github.com/elimuinformatics/vcf2fhir).

## vcfPrepper
Implements the molecular consequence pipeline described on the [Getting Started](https://github.com/FHIR/genomics-operations/wiki/2.-Getting-Started#molecular-consequences) page.
