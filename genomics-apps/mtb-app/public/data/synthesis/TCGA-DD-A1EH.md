# AI Case Synthesis: Higado Sobreviviente (MRN TCGA-DD-A1EH)

> **Demo: synthesis pre-generated for this case** (2026-11-06). Research demonstration, synthetic patient, not for clinical use.

**Contents:** [Reason for visit](#reason-for-visit) · [Case summary](#case-summary) · [Genomic findings](#genomic-findings) · [Cancer drivers](#cancer-drivers) · [Information gaps](#information-gaps) · [Treatment options](#treatment-options) · [Emerging evidence](#emerging-evidence) · [References](#references)

## Reason for visit

Hepatocellular carcinoma with liver and lung recurrence, progressed on sorafenib, now recorded as receiving tazemetostat in an ARID1A trial (NCT05023655) with stable disease (2026-10-15). **That trial is terminated, and tazemetostat was withdrawn worldwide in March 2026 for safety reasons.**

- Should tazemetostat be stopped, and what second-line treatment is supported?
- Do his genomic findings open any other targeted options or trials?
- How should his hepatitis B be managed with the next treatment?

## Case summary

- **Patient:** 24-year-old man. Chronic hepatitis B (diagnosed 2022), on tenofovir since 2026-01. Latex allergy.
- **Diagnosis:** hepatocellular carcinoma of the right lobe on liver biopsy (2025-12-20). AFP 780 ng/mL and Child-Pugh A (score 5) (2025-12-31). Clinical stage IIIA (cT3 cN0 cM0). ECOG 0.
- **Surgery (2026-01-26):** right hepatectomy; pT3 pN0, stage IIIA. Background liver: inactive chronic hepatitis B with fibrosis.
- **Recurrence (2026-05-02):** intrahepatic recurrence in the left lobe and lung metastases. Sorafenib 400 mg twice daily started.
- **Progression (2026-07-14)** on sorafenib.
- **Tumor sequencing (2026-07-22):** ARID1A p.Arg1721\* and p.Val2169fs; JAK1 p.Ala1086Ser; microsatellite stable; tumor mutational burden 1.3 mutations/Mb. ECOG 0 (2026-07-28).
- **Trial treatment:** on study NCT05023655 (tazemetostat for ARID1A-mutated solid tumors) from 2026-08-09; sorafenib stopped at enrollment; tazemetostat 800 mg twice daily recorded as active (order dated 2026-10-15).
- **Current status (2026-10-15):** stable disease. No recent labs, AFP, imaging reports or Child-Pugh score since the recurrence are recorded.

## Genomic findings

Acquired (somatic) tumor changes. The first three are on the 2026-07-22 report; the others come from the tumor's genomic data.

| Gene | Finding | Evidence |
|---|---|---|
| ARID1A | **p.Arg1721\*** (truncating) and **p.Val2169CysfsTer56** (frameshift) | Loss of function of a tumor suppressor; R1721 is a recurrent hotspot (13 tumors). Two truncating changes suggest loss of both copies. CIViC evidence in liver cancer is preclinical only (level D). The ARID1A trial he joined (tazemetostat) is terminated |
| JAK1 | p.Ala1086Ser | Not a known activating change. Report and CIViC: ruxolitinib not supported in hepatocellular carcinoma (level D, preclinical) |
| AXIN1 | **p.Arg647ProfsTer17** (frameshift), not on the report | Loss of function of a tumor suppressor that normally restrains WNT/β-catenin signaling. No ClinVar or CIViC evidence. The one AXIN1-directed trial found (NCT06005974) is terminated |
| IRF2 | **p.Glu191Ter** (truncating), not on the report | Loss of function of a tumor suppressor in interferon signaling. No ClinVar or CIViC evidence |

Also in the genomic data: INHBA p.Arg308Cys, a missense change of unknown significance (not counted). Microsatellite stable with a low tumor mutational burden, so the tissue-agnostic immunotherapy indications (MSI-H, TMB-high) don't apply. No finding carries FDA-labeled or CIViC level A/B therapeutic evidence in hepatocellular carcinoma.

## Cancer drivers

**Bottom line: three brakes are lost (ARID1A, AXIN1, IRF2) and no accelerator was found; none is targetable today, so standard second-line treatment drives the decision.**

| Pathway | Gene | Role | Change | Treatment status |
|---|---|---|---|---|
| Chromatin remodeling | ARID1A | Brake | Two truncating changes | Watch: the EZH2-inhibitor approach (tazemetostat) is withdrawn; another EZH2 inhibitor trial excludes patients who've had one |
| WNT signaling | AXIN1 | Brake | Truncating | Watch: AXIN1-directed trial terminated |
| Interferon signaling | IRF2 | Brake | Truncating | Watch |

ARID1A helps package DNA so genes are switched on and off correctly; its loss was the reason for his tazemetostat trial. AXIN1 normally keeps β-catenin, a WNT-pathway growth signal, in check, so its loss switches that signal on. IRF2 is a brake in the cell's interferon response. None of these is a treatment target outside trials, and the trials aimed at ARID1A and AXIN1 have closed. The data can't show which changes are present in every tumor cell. The Genomics Explorer's Pathways view shows AXIN1 on the WNT pathway diagram.

## Information gaps

- **Tazemetostat status:** the record shows him enrolling (2026-08-09) and receiving tazemetostat after the drug's worldwide withdrawal and the trial's termination (March 2026). Whether he is actually taking it, and for how long, should be confirmed.
- **Current AFP** (needed for ramucirumab eligibility, AFP ≥ 400 ng/mL; also for monitoring). The only value (780 ng/mL) predates resection.
- **Current liver function and Child-Pugh class:** last recorded 2025-12-31, before resection and recurrence. Needed to choose systemic therapy.
- **Extent of disease:** no imaging reports for the 2026-05 recurrence, the 2026-07 progression or the current assessment.
- **Hepatitis B status:** viral load and antigen status aren't recorded; relevant to immunotherapy safety.
- **Genomic data source:** the genomic data's test is dated 2020-01-20, years before the 2026 report, and holds only small variants; variant allele frequencies are recorded as 0, a placeholder.

## Treatment options

**Second line after sorafenib** (immunotherapy-naive, ECOG 0, Child-Pugh A at last assessment, hepatitis B on tenofovir):

**Not recommended: continuing tazemetostat.** The manufacturer withdrew it worldwide on 2026-03-09 and ended all its trials and expanded access after a safety signal (secondary blood cancers); patients on trials were moved to standard care. Its trial for him (NCT05023655) is terminated. Its label warns of secondary malignancies (T-cell lymphoblastic lymphoma, myelodysplastic syndrome, AML, B-cell ALL), so long-term monitoring is advised after exposure.

1. **Immunotherapy** (rank 1, tied; too close to call on the available evidence):
   - **Pembrolizumab:** labeled for HCC secondary to hepatitis B after prior systemic therapy other than a PD-1/PD-L1 regimen, which matches him exactly.
   - **Nivolumab + ipilimumab:** labeled for unresectable or metastatic HCC previously treated with sorafenib.
   - **Durvalumab + tremelimumab:** labeled for unresectable HCC.

   *Why first:* labeled for his setting, and he hasn't had immunotherapy. *Considerations:* continue tenofovir and coordinate hepatitis B management with hepatology (ASCO provisional clinical opinion on HBV in patients receiving cancer therapy); reassess liver function first.
2. **Regorafenib** or **cabozantinib** (rank 2, tied). Each labeled for HCC previously treated with sorafenib. *Why second:* labeled, but another oral multikinase inhibitor after progression on sorafenib.
3. **Ramucirumab** (rank 3). Labeled after sorafenib for HCC with AFP ≥ 400 ng/mL. *Why third:* needs a current AFP; his only value predates surgery.

**Not indicated by his findings:**
- **Atezolizumab + bevacizumab, lenvatinib:** labeled for first-line treatment (atezolizumab + bevacizumab for patients with no prior systemic therapy); he has had sorafenib.
- **Pembrolizumab for MSI-H or TMB-high tumors:** he is microsatellite stable with TMB 1.3 mutations/Mb (pembrolizumab remains an option through its HCC indication, above).
- **Ruxolitinib:** JAK1 p.Ala1086Ser doesn't support sensitivity (preclinical evidence).

**Potential clinical trials** (recruiting on ClinicalTrials.gov as of 2026-11-06; eligibility must be confirmed with each site):

| Trial | Design | Locations |
|---|---|---|
| NCT07042919 | Phase 1b/2: zanzalintinib for advanced HCC in second line and beyond (up to two prior lines, no prior cabozantinib; Child-Pugh A or B7–B8) | Chicago, IL |

Not eligible: NCT04104776 (tulmimetostat, which includes an ARID1A-mutant cohort) excludes prior EZH2 inhibitor treatment and active hepatitis B or C infection; NCT06811116 (SAPHIRE) requires a β-catenin (CTNNB1) mutation and prior immunotherapy; NCT05733598 requires prior PD-1/PD-L1 therapy. Terminated: NCT05023655 (tazemetostat, ARID1A) and NCT06005974 (REC-4881, AXIN1).

## Emerging evidence

Searches covering 2025-10-19 to 2026-11-06 (PubMed and web, FDA announcements, ASCO/ESMO/AACR) found:

- **Tazemetostat withdrawal** (Ipsen, 2026-03-09): withdrawn worldwide in all indications and all trials ended after secondary hematologic malignancies were flagged in the SYMPHONY-1 trial. Directly affects his current treatment (see Treatment options).
- **No new second-line HCC approvals.** Camrelizumab + rivoceranib (first line) received a third FDA complete response letter in July 2026, for manufacturing reasons.

---

## References

- FDA prescribing information (openFDA), retrieved 2026-11-06: Keytruda (pembrolizumab, IV) 2026-07-10; Opdivo (nivolumab) 2026-08-12; Yervoy (ipilimumab) 2026-08-20; Imfinzi (durvalumab) 2026-08-31; Imjudo (tremelimumab) 2026-08-31; Stivarga (regorafenib) 2026-02-09; Cabometyx (cabozantinib) 2025-10-19; Cyramza (ramucirumab) 2026-06-09; Tecentriq (atezolizumab) 2026-05-20; Lenvima (lenvatinib) 2026-06-30; Tazverik (tazemetostat) 2026-04-29, including warnings.
- Ipsen announcement, voluntary worldwide withdrawal of Tazverik (tazemetostat), 2026-03-09 (Targeted Oncology, 2026-03-10; Pharmaceutical Executive, 2026-03-09).
- CIViC (civicdb.org API, retrieved 2026-11-06): JAK1 A1086S 10058 (ruxolitinib, HCC, level D); ARID1A loss 7328 (sorafenib, liver cancer, level D).
- ClinVar (retrieved 2026-11-06): 4530199 (ARID1A p.Arg1721Ter, unclassified). No records for AXIN1 c.1939dup, IRF2 p.Glu191Ter, JAK1 p.Ala1086Ser, INHBA p.Arg308Cys.
- Cancer Hotspots (retrieved 2026-11-06): ARID1A R1721.
- OncoKB Cancer Gene List (downloaded 2026-11-06). TCGA pathways: Sanchez-Vega F, et al. *Cell* 2018;173(2):321-337.e10. Ensembl VEP (GRCh37 REST, 2026-11-06).
- Hepatitis B Virus Screening and Management for Patients With Cancer Prior to Therapy: ASCO Provisional Clinical Opinion Update. *J Clin Oncol* 2020. DOI 10.1200/JCO.20.01757.
- ClinicalTrials.gov, retrieved 2026-11-06. Recruiting: NCT07042919. Not eligible: NCT04104776, NCT06811116, NCT05733598. Terminated: NCT05023655, NCT06005974.
- **Genes searched** (2026-11-06): the tumor's genomic data were searched for small variants in liver cancer genes and in 852 of 858 known cancer genes (six couldn't be located). The genomic test doesn't detect structural or copy-number changes. Trials open to any solid tumor were matched against his alterations on the same date.
