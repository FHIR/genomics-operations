# Case synthesis: Higado Sobreviviente (MRN TCGA-DD-A1EH)

> **Demo: synthesis pre-generated for this case** (2026-10-22). Research demonstration, synthetic patient, not for clinical use.

## Case summary

- **Patient:** 24-year-old Hispanic man with chronic hepatitis B, on tenofovir since January 2026. Latex allergy.
- **Diagnosis (January 2026, age 23):** hepatocellular carcinoma of the right lobe, clinical T3 N0 M0, stage IIIA. AFP 780 ng/mL; Child-Pugh 5 (class A); ECOG 0.
- **Surgery:** right hepatectomy, January 2026 (pathologic T3 N0, stage IIIA). Background liver showed inactive chronic hepatitis with fibrosis.
- **Recurrence (May 2026):** left-lobe liver and lung metastases. Sorafenib started.
- **Clinical progression (July 2026).**
- **Tumor sequencing:** microsatellite stable, TMB 1.3 mutations/Mb. The sequencing report identified candidate trials for ARID1A and AXIN1 alterations.
- **Current treatment:** enrolled in a tazemetostat trial for ARID1A-mutated tumors (NCT05023655). Tazemetostat started and sorafenib stopped in October 2026; disease stable at that visit.
- **Current status:** ECOG 0 (July 2026). Liver function has not been re-recorded since diagnosis.

## Genomic findings

| Gene | Finding | Evidence |
|---|---|---|
| ARID1A | p.Arg1721* (truncating), somatic | No CIViC or ClinVar evidence |
| ARID1A | p.Val2169fs (frameshift), somatic | No CIViC or ClinVar evidence |
| AXIN1 | p.Arg647fs (frameshift), somatic | No CIViC or ClinVar evidence |
| JAK1 | p.Ala1086Ser, somatic | No CIViC or ClinVar evidence |

No knowledge-base therapeutic evidence was found for these variants. The trial matches rest on the gene alterations named in the trials' eligibility, not on CIViC evidence.

## Information gaps

- **Current AFP** (required for ramucirumab eligibility; also useful for monitoring).
- **Current liver function / Child-Pugh class** after hepatectomy and recurrence.
- **Extent of disease** at the July 2026 progression; no imaging reports are recorded.

## Treatment options

**Current plan:** continue on the tazemetostat trial while disease is stable. Tazemetostat is FDA-labeled only for epithelioid sarcoma and follicular lymphoma, so its use in HCC is limited to the trial.

**If disease progresses or the trial ends**, these options are FDA-labeled after sorafenib:
1. **Pembrolizumab.** Labeled for HCC secondary to hepatitis B after prior systemic therapy other than a PD-1/PD-L1-containing regimen. This matches his history.
2. **Regorafenib.** Labeled for HCC previously treated with sorafenib.
3. **Cabozantinib.** Labeled for HCC previously treated with sorafenib.
4. **Ramucirumab.** Labeled after sorafenib for HCC with AFP ≥ 400 ng/mL. His only recorded AFP (780 ng/mL) predates resection, so a current value is needed.
5. **Durvalumab + tremelimumab.** Labeled for unresectable HCC.

*Considerations:*
- **Hepatitis B:** continue antiviral therapy, and coordinate hepatitis B management with hepatology, particularly if an immune checkpoint inhibitor is used. HBV management in patients receiving cancer therapy is addressed by an ASCO provisional clinical opinion.
- **Liver function:** current function should be reassessed before choosing systemic therapy.

**Not indicated by his findings:**
- **Nivolumab + ipilimumab:** labeled for HCC as first-line treatment only.
- **Pembrolizumab for TMB-high or MSI-high tumors:** requires TMB ≥ 10 mutations/Mb or MSI-H/dMMR. He is MSS with TMB 1.3.

**Potential clinical trials** (ClinicalTrials.gov, checked 2026-10-22):

| Trial | Title | Status |
|---|---|---|
| NCT05023655 | Phase II Study of Tazemetostat in Solid Tumors Harboring an ARID1A Mutation | Terminated on ClinicalTrials.gov; enrollment in this scenario is illustrative |
| NCT06005974 | A Study of REC-4881 in Participants With Cancers Which Have an AXIN1 or APC Mutation | Terminated |

---

## References

- FDA prescribing information (openFDA), retrieved 2026-10-22: Keytruda (pembrolizumab, IV), label 2026-07-10; Stivarga (regorafenib) 2026-02-09; Cabometyx (cabozantinib) 2025-10-19; Cyramza (ramucirumab) 2026-06-09; Imjudo (tremelimumab) 2026-08-31; Imfinzi (durvalumab) 2026-08-31; Opdivo (nivolumab) 2026-08-12; Yervoy (ipilimumab) 2026-08-20; Tazverik (tazemetostat) 2026-04-29.
- ClinicalTrials.gov, retrieved 2026-10-22: NCT05023655, NCT06005974.
- Hepatitis B Virus Screening and Management for Patients With Cancer Prior to Therapy: ASCO Provisional Clinical Opinion Update. *J Clin Oncol* 2020. DOI 10.1200/JCO.20.01757.
