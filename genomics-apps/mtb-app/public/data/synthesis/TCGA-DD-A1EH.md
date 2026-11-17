# AI Case Synthesis: Higado Sobreviviente (MRN TCGA-DD-A1EH)

> **Demo: synthesis pre-generated for this case** (2026-11-16). Research demonstration, synthetic patient, not for clinical use.

**Contents:** [Reason for visit](#reason-for-visit) · [Case summary](#case-summary) · [Genomic findings](#genomic-findings) · [Cancer drivers](#cancer-drivers) · [Information gaps](#information-gaps) · [Treatment options](#treatment-options) · [Emerging evidence](#emerging-evidence) · [References](#references)

## Reason for visit

Hepatitis B-related hepatocellular carcinoma with liver and lung recurrence that progressed on sorafenib (2026-07-14). The record shows enrollment in trial NCT05023655 (2026-08-09) and tazemetostat started 2026-10-15, but ClinicalTrials.gov lists that trial as terminated on 2026-03-06, and the FDA has since warned against tazemetostat and withdrawn its approval.

- Should tazemetostat be continued?
- What next-line therapy is supported after sorafenib?
- Do the tumor's ARID1A, AXIN1 and IRF2 changes open any targeted option or trial?

## Case summary

- **Patient:** 24-year-old man (born 2002-10-14). Chronic hepatitis B since 2022, on tenofovir disoproxil fumarate 300 mg daily since 2026-01-06. Latex allergy (high criticality).
- **Diagnosis:** hepatocellular carcinoma (HCC) of the right liver lobe.
  - Liver biopsy 2025-12-20. Alpha-fetoprotein (AFP) 780 ng/mL and Child-Pugh A (score 5) on 2025-12-31.
  - CT 2026-01-06: two LI-RADS 5 masses in the right lobe (7.8 cm, segment 7; 2.1 cm, segment 6); no portal or hepatic vein invasion; features of chronic liver disease; no extrahepatic disease. Clinical stage IIIA (cT3 cN0 cM0).
- **Treatment:**
  - Right hepatectomy 2026-01-26 (curative intent): two HCCs (7.8 and 2.1 cm), margins negative; background chronic hepatitis, inactive, with fibrosis. Pathologic stage IIIA (pT3 pN0). CT 2026-03-11: no residual or recurrent tumor.
  - Sorafenib 400 mg twice daily from 2026-05-02 (palliative), stopped 2026-10-15 for absent response.
  - Tazemetostat 800 mg twice daily from 2026-10-15, on trial NCT05023655 (on study since 2026-08-09).
- **Recurrence and progression:**
  - CT 2026-05-02: new 2.6 cm lesion in the left lobe (segment 3) and multiple new bilateral lung nodules (largest 1.4 cm), consistent with recurrence and lung metastases.
  - CT 2026-07-13: lung metastases larger and more numerous (largest 2.2 cm); left-lobe recurrence 3.9 cm. MRI 2026-07-14: left-lobe recurrence 4.0 cm, new 1.2 cm nodule in segment 2, no portal vein tumor thrombus, no ascites. Progression recorded 2026-07-14.
- **Tumor testing** (dated 2026-07-22, reported 2026-07-28): ARID1A p.Arg1721\* and p.Val2169fs; JAK1 p.Ala1086Ser; microsatellite stable; tumor mutational burden (TMB) 1.3 mutations/Mb.
- **Current status:** ECOG 0 (2026-07-28). Disease recorded as stable on 2026-10-15, with no imaging report linked to that assessment.

## Genomic findings

All findings are acquired (somatic) changes in the tumor. No germline (inherited) testing is in the record.

| Gene | Finding | Evidence |
|---|---|---|
| ARID1A | **p.Arg1721\*** (c.5161C>T), stop-gained | Reported, and present in the genomic data. A recurrent mutation site in tumors (Cancer Hotspots). ClinVar: one somatic record, clinical impact "Tier II – potential"; no oncogenicity classification. CIViC: no evidence for this variant, and none at level A or B for a therapy in HCC (ARID1A items are preclinical, functional, or in other cancers such as ovarian, kidney, stomach and lymphoma) |
| ARID1A | **p.Val2169Cysfs\*56** (c.6504_6505insT), frameshift | Reported, and present in the genomic data. Not in ClinVar or CIViC |
| AXIN1 | **p.Arg647Profs\*17** (c.1939dup), frameshift | **Not on the genetic report**; found in the genomic data. Not in ClinVar, CIViC or Cancer Hotspots |
| IRF2 | **p.Glu191\*** (c.571G>T), stop-gained | **Not on the genetic report**; found in the genomic data. Not in ClinVar, CIViC or Cancer Hotspots |
| JAK1 | p.Ala1086Ser (c.3256G>T), missense | Reported, and present in the genomic data. CIViC (preclinical, level D): does not activate JAK1 and does not support ruxolitinib sensitivity in HCC. Not a recurrent mutation site; not in ClinVar |
| INHBA | p.Arg308Cys (c.922C>T), missense | Not on the genetic report; found in the genomic data. Effect unknown: not in ClinVar, CIViC or Cancer Hotspots. No therapeutic or trial relevance |

Both ARID1A changes truncate the protein. Whether they hit the same copy of the gene or both copies can't be told from this data. The genomic data hold three changes the report doesn't list: truncating changes in AXIN1 and IRF2, and an INHBA missense change of unknown effect; the lab should confirm or explain them. No other small variant was found across all known cancer genes. Within the regions the test studied, that rules out the common β-catenin (CTNNB1) hotspot mutations, TP53 mutations and BRAF V600E. The TERT promoter hotspots lie outside the studied regions and weren't tested. The test detects only small variants, so amplifications (e.g. FGF19) and fusions (e.g. NTRK, RET, FGFR2) weren't tested.

## Cancer drivers

**Bottom line: the tumor is driven by loss of three brakes, ARID1A (a chromatin regulator, hit twice), AXIN1 (a brake on the WNT/β-catenin growth pathway) and IRF2 (an interferon-response regulator); none is targetable with an approved drug today.**

| Pathway | Gene | Role | Change | Treatment status |
|---|---|---|---|---|
| Chromatin remodeling (SWI/SNF) | ARID1A | Brake | Two truncating changes (stop-gained, frameshift) | Trial only (see Treatment options) |
| WNT/β-catenin | AXIN1 | Brake | Truncating frameshift | Trial only (needs confirmation) |
| Interferon signaling | IRF2 | Brake | Truncating (stop-gained) | Watch |

ARID1A helps package DNA so that the right genes are switched on and off; with it broken, the cell's growth controls loosen. AXIN1 holds β-catenin in check, so losing it lets the WNT growth pathway run, typically more weakly than the β-catenin mutations common in liver cancer. IRF2 is a transcription factor in the interferon response that also restrains cell growth; its loss is a recognized, though uncommon, event in liver cancer, with no drug or trial linked to it. All three are typical of HCC. The JAK1 and INHBA changes aren't drivers: they're missense changes of unknown effect, and for JAK1 laboratory evidence says it doesn't switch the gene on. Worth confirming: that the AXIN1 and IRF2 changes are real (they aren't on the report), and whether ARID1A protein is lost (immunohistochemistry). Clonality can't be judged, because the variant allele fractions in the data are placeholders. The app's Pathways view shows the WNT pathway diagram.

## Information gaps

- **Which treatment the patient is on now.** The trial he is listed on closed on 2026-03-06, before his recorded enrollment (2026-08-09); tazemetostat start (2026-10-15) and sorafenib stop (2026-10-15) come three months after progression. The record needs to say what he is taking, and under what protocol or access program.
- **Current disease status.** The last imaging is from 2026-07-14; the 2026-10-15 "stable" assessment has no imaging linked. Current measurable disease is needed for any trial.
- **Current liver function and AFP.** Child-Pugh (A5) and AFP (780 ng/mL) date from 2025-12-31, before resection. The trials behind each option below enrolled Child-Pugh A, and the listed trials require Child-Pugh A or, in several, no worse than B7; ramucirumab needs a current AFP ≥400 ng/mL.
- **Hepatitis B viral load.** Not recorded. Pembrolizumab's supporting trial required controlled HBV (viral load <2000 IU/mL), and several trials require <500 IU/mL or undetectable.
- **Recent labs and ECOG.** No blood counts, bilirubin, albumin or INR, and ECOG (0) is from 2026-07-28.
- **Tumor sample for sequencing.** The report doesn't say which sample was tested on 2026-07-22; no biopsy is recorded on that date. The test date in the genomic data (2020-01-20) also doesn't match the report.
- **AXIN1, IRF2 and INHBA changes** aren't on the genetic report; the lab should confirm or explain them. The AXIN1 change matters for one trial's eligibility.
- **TERT promoter** (mutated in many liver cancers) lies outside the regions the test studied, so it wasn't tested.
- **Copy-number changes and rearrangements** weren't tested: e.g. FGF19 amplification (a trial target in HCC), MET or ERBB2 amplification, NTRK, RET or FGFR2 fusions. RNA-based fusion testing or a copy-number-capable panel would answer this.
- **Protein tests used by trials:** glypican-3 (GPC3) for two trials; FGF19 for one; ARID1A protein loss for ATR-inhibitor studies; PD-L1.
- **Germline testing:** none recorded. He was diagnosed at 23; chronic hepatitis B is a known cause, but whether inherited risk was considered isn't recorded.

## Treatment options

### Not recommended: continuing tazemetostat

In May 2026 the FDA alerted clinicians to an increased rate of new blood cancers with tazemetostat (18 of 318 patients, 5.7%, in its confirmatory trial, mostly myelodysplastic syndrome and acute myeloid leukemia, some fatal), stated that its risks outweigh its benefits, and reported that the independent monitoring committee recommended patients stop it; the approval was withdrawn (Federal Register, 2026-06-22). The sponsor discontinued all tazemetostat trials, and NCT05023655 is listed as terminated. Its label (Tazverik, 2026-04-29) never covered liver cancer (epithelioid sarcoma and follicular lymphoma only), and CIViC has no evidence for an EZH2 inhibitor in ARID1A-mutant HCC. **Stop tazemetostat** and move to one of the options below; because new blood cancers have appeared after stopping the drug, blood-count monitoring is reasonable.

### Next-line therapy after sorafenib (ranked)

**1. Pembrolizumab (Keytruda IV, or Keytruda Qlex subcutaneous)**
- *Label:* HCC secondary to hepatitis B in patients who have had prior systemic therapy other than a PD-1/PD-L1 regimen (Keytruda and Keytruda Qlex labels 2026-10-16). He fits: hepatitis B-related HCC, prior sorafenib, no prior immunotherapy.
- *Evidence:* KEYNOTE-394 (phase 3, placebo-controlled, after sorafenib or oxaliplatin chemotherapy): overall survival 13.9 vs 13.0 months (hazard ratio 0.78); response rate 11% vs 1.6%, median response duration 23.9 months.
- *Considerations:* the trial required Child-Pugh A and controlled hepatitis B (viral load <2000 IU/mL), so a current viral load and Child-Pugh are needed; continue tenofovir. ECOG 0 and no autoimmune disease recorded. Microsatellite-stable status and low TMB don't exclude him from this disease-specific label. Exploratory first-line data suggest a smaller advantage of PD-1 therapy over sorafenib when the WNT pathway is altered (see Emerging evidence); not a reason to withhold it.
- *Why first:* a label written for exactly his setting, a survival benefit in a phase 3 trial, a mechanism he hasn't been exposed to, and durable responses in the minority who respond.

**2. Regorafenib (Stivarga) or cabozantinib (Cabometyx), ranked together**
- *Labels:* both for HCC previously treated with sorafenib (Stivarga 2026-02-09; Cabometyx 2025-10-19).
- *Evidence:* RESORCE (phase 3, after progression on sorafenib): overall survival 10.6 vs 7.8 months (hazard ratio 0.63). CELESTIAL (phase 3, after sorafenib, Child-Pugh A): 10.2 vs 8.0 months (hazard ratio 0.76).
- *Considerations:* RESORCE required tolerance of sorafenib at ≥400 mg/day; he took 400 mg twice daily for about five months, with no recorded intolerance. Both trials enrolled Child-Pugh A; Stivarga's label warns of liver toxicity, bleeding and high blood pressure. Either is a reasonable sequel if pembrolizumab fails or can't be given.
- *Why second:* same strength of evidence as pembrolizumab, but both are kinase inhibitors that overlap with the sorafenib he progressed on.

**3. Ramucirumab (Cyramza), only if AFP ≥400 ng/mL**
- *Label:* HCC with AFP ≥400 ng/mL after sorafenib (Cyramza label 2026-06-09). REACH-2 (phase 3): overall survival 8.5 vs 7.3 months (hazard ratio 0.71).
- *Considerations:* his only AFP (780 ng/mL) predates resection; a current value is needed.
- *Why third:* a phase 3 survival benefit, but eligibility depends on a value not yet known, and the benefit is the smallest of the randomized options.

**4. Nivolumab + ipilimumab (Opdivo + Yervoy)**
- *Label:* unresectable or metastatic HCC previously treated with sorafenib (Opdivo 2026-08-12; Yervoy 2026-08-20).
- *Evidence:* CheckMate-040 cohort 4 (single-arm, 49 patients who progressed on or were intolerant to sorafenib, Child-Pugh A; 57% had active hepatitis B): response rate 33%, with 56% of responses lasting at least 12 months. No randomized survival data in this setting.
- *Considerations:* needs Child-Pugh A and no autoimmune disease; combination immunotherapy carries more immune-related toxicity than pembrolizumab alone. It's an alternative immunotherapy to pembrolizumab, not a sequel to it: the trial didn't enroll patients after PD-1 therapy.
- *Why fourth:* labeled for his setting, but supported only by response rate in a small single-arm cohort, behind the options with phase 3 survival benefit.

Clinical trials (below) are a reasonable alternative to any of these and can be considered alongside them; none ranks above the labeled options on current evidence.

The co-occurring-variant evidence in CIViC was checked for his altered genes (ARID1A, AXIN1, IRF2, JAK1, INHBA); no combination evidence applies.

### Not indicated by the patient's findings

- **Atezolizumab + bevacizumab, lenvatinib:** labeled for untreated or first-line HCC only (Tecentriq 2026-10-13; Lenvima 2026-10-14). He has had sorafenib. (Nivolumab + ipilimumab is labeled both first-line and after sorafenib; see option 4.)
- **Durvalumab + tremelimumab:** the label says "unresectable HCC" without a line of therapy (Imfinzi and Imjudo 2026-08-31), but its supporting trial (HIMALAYA) enrolled only patients with no prior systemic therapy, so evidence after sorafenib is lacking.
- **Ruxolitinib for JAK1 p.Ala1086Ser:** CIViC preclinical evidence (level D) says this variant doesn't activate JAK1 and doesn't support ruxolitinib sensitivity in HCC; the lab report says the same.
- **Tissue-agnostic immunotherapy (MSI-H, TMB-H):** microsatellite stable, TMB 1.3 mutations/Mb.
- **Tissue-agnostic targeted drugs (NTRK, RET, BRAF V600E, HER2):** no BRAF V600E or other qualifying small variant was found; fusions and HER2 weren't tested (see Information gaps).
- **Sorafenib rechallenge:** progressed on it.

### Potential clinical trials

Searched on ClinicalTrials.gov 2026-11-16 (recruiting, US sites): trials for ARID1A, AXIN1, IRF2 and WNT-pathway alterations, trials open to any solid tumor that were matched against his alterations, and HCC trials after prior systemic therapy. Eligibility must be confirmed with each site; all would need a washout from tazemetostat.

| Trial | Design | Locations |
|---|---|---|
| NCT05440708 | TTI-101 (STAT3 inhibitor) alone (Cohort A: progression on up to 3 prior lines), phase 1/2. HBV allowed if viral load <500 IU/mL on antivirals | 15 US sites incl. Tampa FL, Baltimore MD, Houston and Dallas TX, New York NY |
| NCT05797805 | Tegavivint (TBL1 inhibitor, acts on β-catenin signaling) alone or with cabozantinib or lenvatinib, phase 1/2, after ≥1 prior line. HBV allowed if viral load <500 IU/mL | Charlotte NC, Miami FL, Chicago IL, Houston and Dallas TX, Seattle WA, Duarte CA |
| NCT06600321 | ALN-BCAT (RNA therapy against β-catenin) alone or with pembrolizumab, phase 1. Requires a WNT-pathway activating mutation: the site must confirm whether his AXIN1 frameshift qualifies | 15 US sites incl. Atlanta GA, Jacksonville FL, Richmond VA |
| NCT06858813 | ABBV-324 (Part 1: after failure of ≥1 prior systemic therapy), phase 1 | 10 US sites incl. Nashville TN, Philadelphia PA, Chicago IL |
| NCT07042919 | Zanzalintinib (kinase inhibitor) after 1–2 prior lines (no prior cabozantinib), phase 1/2 | Chicago IL |
| NCT07174570 | Celecoxib + durvalumab + tremelimumab, phase 2 (no immunotherapy in prior 6 months) | Atlanta GA |
| NCT04864054 | ECT204 (GPC3-targeted T cells) after ≥ first-line therapy, phase 1/2. Needs GPC3-positive tumor by IHC; treated HBV allowed | Dallas TX, Seattle WA, Portland OR, Duarte CA, Bronx NY |
| NCT06764316 | BAY 3547926 (GPC3-targeted), phase 1. Needs GPC3-positive tumor by IHC | Dothan AL, Duarte CA |
| NCT06915753 | TYRA-430 (FGF/FGFR pathway drug), phase 1. HCC after standard care (Part A, no genomic requirement; Part B HCC cohort needs archival tissue for central FGF19 IHC). Child-Pugh A | Los Angeles, San Francisco and Stanford CA, Westwood KS, Baltimore MD, Boston MA, Detroit MI, New York NY |
| NCT05490472 | JAB-2485 (Aurora A inhibitor), phase 1/2, expansion includes ARID1A-mutant solid tumors. Excludes "active" hepatitis B: needs the site's view of treated, suppressed HBV | Detroit MI, St Louis MO, Dallas TX, Salt Lake City UT |

Not eligible: NCT04104776 (tulmimetostat, ARID1A-mutant solid tumors) excludes prior EZH2 inhibitors, which tazemetostat is; NCT06726161 (GPC3 radioligand) excludes prior participation in an interventional study; NCT06811116 (sapanisertib + cabozantinib) requires prior checkpoint-inhibitor therapy unless ineligible, and its phase 2 part a β-catenin (CTNNB1) mutation, which wasn't found; NCT05733598, NCT07219459 and NCT05168163 require prior immunotherapy.

## Emerging evidence

Searches covered 2025-10-19 to 2026-11-16 (PubMed, FDA announcements, web).

- **Tazemetostat withdrawn for safety** (FDA safety alert, May 2026; Ipsen announcement 2026-03-09; approval withdrawn, Federal Register 2026-06-22). Confirmatory-trial finding of new blood cancers in 5.7% of treated patients. Bears directly on the decision to stop it (see Treatment options).
- **Tazemetostat in SWI/SNF-altered solid tumors** (phase 2 basket, *Nat Commun* 2026): responses of 0–9% in rhabdoid, INI1-negative, synovial sarcoma, renal medullary and chordoma cohorts; no ARID1A or liver cancer cohort. Gives no support for EZH2 inhibition in his tumor.
- **ATR inhibitor ceralasertib in ARID1A-deficient tumors** (phase 2, *Clin Cancer Res* 2026): response rate 14% overall, all responses in endometrial or ovarian clear-cell cancers; eligibility was ARID1A protein loss by IHC. No liver cancer responses reported; a possible future direction only.
- **Genomically matched therapy in refractory liver cancer** (French genomic medicine initiative, *J Hepatol* 2026): disease control in 23% of HCC patients given matched therapy, and only for alterations with stronger evidence tiers. Tempering for options built on ARID1A or AXIN1 alone.
- **WNT pathway and PD-1 benefit** (CheckMate 459 exploratory biomarker analysis, *JHEP Rep* 2026): sorafenib did better in tumors with WNT/β-catenin alterations, so nivolumab's survival advantage over sorafenib was confined to tumors without them. Exploratory and first-line; it bears on how immunotherapy compares with kinase inhibitors given his AXIN1 change, but isn't grounds to withhold pembrolizumab.
- **AXIN1 mutations in HCC** (*iScience* 2026, laboratory study): AXIN1 truncations activate WNT signaling, but generally less than β-catenin mutations, and favor YAP/TAZ signaling. Relevant to whether WNT-directed trials apply.
- No new FDA approval for HCC after prior systemic therapy was found.

---

## References

**FDA prescribing information** (openFDA, retrieved 2026-11-16)
- Keytruda (pembrolizumab), label 2026-10-16; Keytruda Qlex (pembrolizumab and berahyaluronidase alfa), label 2026-10-16
- Stivarga (regorafenib), label 2026-02-09
- Cabometyx (cabozantinib), label 2025-10-19
- Cyramza (ramucirumab), label 2026-06-09
- Tecentriq (atezolizumab), label 2026-10-13
- Opdivo (nivolumab), label 2026-08-12; Yervoy (ipilimumab), label 2026-08-20
- Imfinzi (durvalumab) and Imjudo (tremelimumab), labels 2026-08-31
- Lenvima (lenvatinib), label 2026-10-14
- Tazverik (tazemetostat), label 2026-04-29

**FDA and sponsor announcements** (retrieved 2026-11-16)
- FDA. FDA Alerts Health Care Providers and Patients about Increased Risk of New Blood Cancers with Tazverik (tazemetostat) Use; Sponsor to Voluntarily Withdraw Product from Market. May 2026 (page current as of 2026-05-11).
- Ipsen. Update: Ipsen voluntarily withdraws Tazverik (tazemetostat) in follicular lymphoma and epithelioid sarcoma. 2026-03-09.
- Federal Register. Withdrawal of approval of NDA 211723 (Tazverik). Document 2026-12367, 2026-06-22.

**CIViC** (retrieved 2026-11-16)
- JAK1 p.Ala1086Ser: [civicdb.org/evidence/10055](https://civicdb.org/evidence/10055) (level D, does not support gain of function), [civicdb.org/evidence/10058](https://civicdb.org/evidence/10058) (level D, HCC, does not support ruxolitinib sensitivity)
- ARID1A: 14 evidence items reviewed; none at level A or B for a therapy in HCC

**Co-occurring variants:** CIViC co-occurring-variant profiles (retrieved 2026-10-02): none involve ARID1A, AXIN1, IRF2, JAK1 or INHBA.

**CIViC** for AXIN1, IRF2 and INHBA (retrieved 2026-11-16): no evidence items.

**ClinVar** (retrieved 2026-11-16): ClinVar records for ARID1A c.5161C>T (somatic clinical impact Tier II – potential); no records for ARID1A c.6504_6505insT, AXIN1 c.1939dup, IRF2 c.571G>T, INHBA c.922C>T or JAK1 c.3256G>T.

**Cancer genes and hotspots:** OncoKB Cancer Gene List (downloaded 2026-11-06); Cancer Hotspots (retrieved 2026-11-16). Variant effects computed with Ensembl VEP (GRCh37, 2026-11-16), matching the genomic data's own annotations.

**Pathways:** Sanchez-Vega F, et al. Oncogenic Signaling Pathways in The Cancer Genome Atlas. *Cell*. 2018;173(2):321-337.e10. doi:10.1016/j.cell.2018.03.035

**Literature**
- Gounder M, et al. Tazemetostat, an EZH2 inhibitor, in solid tumors harboring SWI/SNF alterations: a phase II basket study. *Nat Commun*. 2026;17:7120. doi:10.1038/s41467-026-69708-2
- Zhu X, et al. Efficacy of the ATR Inhibitor Ceralasertib in Patients with ARID1A-Deficient Gynecologic and Other Solid Tumor Malignancies. *Clin Cancer Res*. 2026;32(9):1632-1640. doi:10.1158/1078-0432.CCR-25-4043
- Campani C, et al. Molecularly guided therapies for advanced primary liver cancers refractory to systemic treatment: Results from the 2025 French genomic medicine initiative. *J Hepatol*. 2026;85(2):289-300. doi:10.1016/j.jhep.2026.03.019
- Biomarkers associated with efficacy of first-line nivolumab vs. sorafenib in advanced hepatocellular carcinoma: Results from the CheckMate 459 phase III trial. *JHEP Rep*. 2026;8(12):102038. doi:10.1016/j.jhepr.2026.102038
- Hepatocellular carcinoma-linked AXIN1 mutations drive low Wnt/β-catenin activity enabling niche-independent growth and YAP/TAZ signaling. *iScience*. 2026;29(1):114501. doi:10.1016/j.isci.2025.114501

**ClinicalTrials.gov** (searched 2026-11-16; status as of that date): NCT05023655 (terminated 2026-03-06); recruiting: NCT05440708, NCT05797805, NCT06600321, NCT06858813, NCT07042919, NCT07174570, NCT04864054, NCT06764316, NCT06915753, NCT05490472; not eligible: NCT04104776, NCT06726161, NCT06811116, NCT05733598, NCT07219459, NCT05168163.

**Genes searched** (2026-11-16): small variants in the liver cancer genes, the genes behind tumor-agnostic approvals, and all 858 known cancer genes (869 genes in all), plus the TERT promoter region; the regions the test studied were checked for the liver cancer genes, the tumor-agnostic genes and every gene stated as unaltered. The test detects only small variants, so copy-number changes and rearrangements couldn't be searched, and the TERT promoter hotspots were outside the studied regions. Three variants in cancer genes beyond those reported were found and are listed above; none were left out.
