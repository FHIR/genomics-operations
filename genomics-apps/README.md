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
* find-subject-intersecting-variants
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
