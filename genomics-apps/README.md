# Genomics Applications
Here are some simple genomics applications built using the FHIR Genomics Operations:

## Simple App
https://simplegenomics.streamlit.app/

Enter patient and genomic range. Raw output from find-subject-variants is shown in left column, and tabular output is shown in right column.

Uses:
* find-subject-variants

## Get Variants
https://getvariants.streamlit.app/

Enter patient and gene, and get back all overlapping simple variants, structural variants, and genotypes.

Uses:
* get-feature-coordinates
* find-subject-variants
* find-subject-structural-intersecting-variants
* find-subject-haplotypes

## Genetic Screening
https://geneticscreening.streamlit.app/

Select a condition from dropdown. The population is screened for this condition, and each associated variant in each patient is returned.

Uses:
* find-population-dx-implications
* find-subject-dx-implications

## PGx Screening
https://pgxscreening.streamlit.app/

Select patient from dropdown. List of potential drug-gene interactions are shown in right column. The patient's med list is shown in left column, where those meds with PGx interactions are flagged.

Uses:
* find-subject-haplotypes
* find-subject-tx-implications

## Mendelian Screening
https://mendelianscreening.streamlit.app/

Enter patient and gene(s), and get back all variants, annotated with diagnostic implications, predicted molecular consequences, and population
allele frequencies. Patient's problem list is flagged with conditions potentially having an underlying genetic basis.

Uses:
* FHIR search (to retrieve conditions, medications, allergies)
* find-subject-variants
* find-subject-molecular-consequences
* find-subject-dx-implications

## Get Polygenic Score
https://polygenicscore.streamlit.app/ (and video: https://vimeo.com/930060292)

This proof-of-concept app demonstrates the computation of a polygenic score, using FHIR Genomics Operations to access a person's entire genome.
Select patient and polygenic model in the sidebar, optionally adjust standard deviation threshold, and click 'Run'.

Uses:
* find-subject-variants
* find-subject-specific-variants
* find-study-metadata
