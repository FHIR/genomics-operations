# Search log: TCGA-DD-A1EH

Lookups made while writing `TCGA-DD-A1EH.md` on 2026-11-06, following `README.md`. Genomics server: https://fhir-gen-ops.herokuapp.com (subject TCGA-DD-A1EH).

## Clinical record

| Query | Result |
|---|---|
| `clinical-data/TCGA-DD-A1EH/TCGA-DD-A1EH-transaction-bundle.json` | Read in full. No referral; Encounters have no reason. Reason for visit inferred from trial enrollment (ResearchSubject on-study 2026-08-09; tazemetostat MedicationRequest active, authored 2026-10-15) and stable disease 2026-10-15. One genetic report 2026-07-22 (ARID1A p.Arg1721\*, p.Val2169fs; JAK1 p.Ala1086Ser; MSS; TMB 1.3; ruxolitinib level D) |

## Genome

| Query | Why | Result |
|---|---|---|
| `$find-study-metadata` | Tests | TCGA-DD-A1EH-T1 (2020-01-20, specimen TCGA-DD-A1EH-Sp1, GRCh37; SNV, MNV, delins only) |
| `$find-subject-variants`, GRCh37, 852 genes in batches of 20; structural search skipped (not detected by the test) | Minimum genes (liver cancer gene list, reported genes, always-search list) and cancer drivers | No errors. 6 variants in 5 cancer genes, all somatic; VAF recorded as 0 |

## Molecular consequences (Ensembl VEP, GRCh37, 2026-11-06)

ARID1A p.Arg1721Ter (stop gained); ARID1A p.Val2169CysfsTer56 (frameshift); JAK1 p.Ala1086Ser (missense, splice region); IRF2 p.Glu191Ter (stop gained); INHBA p.Arg308Cys (missense); AXIN1 p.Arg647ProfsTer17 (frameshift).

## Driver rule (README §2.5)

| Finding | Basis | Driver? |
|---|---|---|
| ARID1A p.Arg1721Ter | TSG, truncating; hotspot R1721 (13 tumors) | Yes |
| ARID1A p.Val2169CysfsTer56 | TSG, truncating (second hit) | Yes |
| AXIN1 p.Arg647ProfsTer17 | TSG, truncating | Yes |
| IRF2 p.Glu191Ter | TSG, truncating | Yes |
| JAK1 p.Ala1086Ser | ONCOGENE_AND_TSG; no hotspot or ClinVar record | No |
| INHBA p.Arg308Cys | ONCOGENE_AND_TSG; no hotspot or ClinVar record | No |

## Knowledge bases

| Query | Result |
|---|---|
| ClinVar | 4530199 ARID1A p.Arg1721Ter (no classification). No records: AXIN1 c.1939dup, IRF2 p.Glu191Ter, INHBA p.Arg308Cys, JAK1 p.Ala1086Ser |
| CIViC | JAK1 A1086S: EID10058 (ruxolitinib, HCC, level D, does not support). ARID1A loss: EID7328 (sorafenib, liver cancer, level D). AXIN1: nothing in liver cancer |
| TCGA pathways | WNT (AXIN1) |
| Crossref, DOI 10.1200/JCO.20.01757 | Hepatitis B Virus Screening and Management for Patients With Cancer Prior to Therapy: ASCO PCO Update, *J Clin Oncol* 2020 |

## Labels (openFDA, retrieved 2026-11-06)

Keytruda IV (2026-07-10): HCC secondary to hepatitis B after prior non-PD-1/PD-L1 systemic therapy. Opdivo (2026-08-12) + Yervoy (2026-08-20): first-line, and after sorafenib. Imfinzi + Imjudo (2026-08-31): unresectable HCC. Stivarga (2026-02-09), Cabometyx (2025-10-19): after sorafenib. Cyramza (2026-06-09): after sorafenib, AFP ≥ 400. Tecentriq (2026-05-20): with bevacizumab, no prior systemic therapy. Lenvima (2026-06-30): first line. Tazverik (2026-04-29): epithelioid sarcoma, follicular lymphoma; warnings: secondary malignancies (T-LBL, MDS, AML, B-ALL).

## Trials (ClinicalTrials.gov API v2, 2026-11-06)

| Search | Result |
|---|---|
| NCT05023655 status | **TERMINATED**: "Ipsen has discontinued all active tazemetostat clinical trials and expanded access programs due to the emergence of a clinically unfavorable benefit risk profile ..." (last update 2026-03-17) |
| NCT06005974 status | TERMINATED (sponsor decision, not safety) |
| HCC, second line / previously treated, recruiting, US | NCT07042919 zanzalintinib (eligible). NCT05733598 RP2 (requires prior PD-1/PD-L1; not eligible). NCT06811116 SAPHIRE (β-catenin mutation and prior ICI; not eligible) |
| ARID1A, recruiting, US | NCT04104776 tulmimetostat cohort M1 (excludes prior EZH2 inhibitor and active hepatitis B/C); others disease-specific (ovarian, endometrial, gastroesophageal) |
| AXIN1, recruiting, US | None |
| Tumor-agnostic checklist (`tumor-agnostic-trials.md`) for ARID1A, AXIN1, IRF2, JAK1 | ARID1A: NCT04104776 only (not eligible) |

## Recent evidence

Search window: 2025-10-19 (oldest label used, Cabometyx) to 2026-11-06.

| Search | Result |
|---|---|
| Ipsen tazemetostat withdrawal | 2026-03-09: voluntary worldwide withdrawal; all trials and expanded access ended; SYMPHONY-1 secondary hematologic malignancies; enrolled patients moved off drug |
| FDA approvals HCC 2025–2026, second line | None. Camrelizumab + rivoceranib (first line): third CRL July 2026 (manufacturing) |
