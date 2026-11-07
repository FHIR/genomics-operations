# Creating a clinical case for the MTB app

How to build a new patient case for the MTB app so that the clinical record, the genetic reports and the genomic data on the server tell one consistent story. Gathered 2026-11-06 from the four existing cases (HCC1143, L2345, MRN1234, TCGA-DD-A1EH) and the problems found in them while writing their case syntheses.

A case has three parts, and they must agree with each other:

1. **The clinical record:** an mCODE transaction bundle in this folder, loaded to the clinical FHIR server (the Case Summary tab).
2. **The genomic data:** the patient's variants on the genomics operations server, loaded from VCFs (the Genomics Explorer tab).
3. **The case synthesis** (optional): written afterward from parts 1 and 2, following `../public/data/synthesis/README.md`.

## 1. Identity: the MRN

The MRN is the key that joins everything. Use the same string everywhere.

| Where | How the MRN appears |
|---|---|
| Patient resource | `identifier` with `system` `http://elimu.io/fhir/mrn`, `type` v2-0203 `MR` ("Medical Record Number"), `use` `usual`, `value` = the MRN. The app finds patients only by this system (`clinicalService.ts`, `MRN_IDENTIFIER_SYSTEM`) |
| Patient id | `mtb-<mrn in lower case>` (e.g. `mtb-l2345`) |
| Bundle id | `mtb-<mrn>-transaction` |
| Other resource ids | `mtb-<mrn>-<description>` (e.g. `mtb-l2345-genomic-variant-egfr-t790m`) |
| Genomics server | `patient_id` in `vcfData.csv` and `patientID` in the tests list = the MRN exactly; the operations' `subject` parameter is the MRN |
| Folder and files | `clinical-data/<MRN>/` with `<MRN>-transaction-bundle.json`, `build_<mrn>.py` (preferred: generate the bundle from a script) and `<MRN>-scenario.md` |
| Synthesis | `public/data/synthesis/<MRN>.md` |

- **One patient per MRN** on the clinical server. If two Patients carry the same MRN, the app warns and uses the first.
- The app's MRN box lists every Patient on the clinical server with an identifier in this system, so a test Patient with a stray MRN shows up there.
- TCGA-DD-A1EH uses `use: official` and the display "Medical record number". The other three use `usual` and "Medical Record Number". Use the latter for new cases.

## 2. The clinical bundle

Conventions:

- **Transaction bundle**, every entry `PUT` to `<Type>/<id>`, with relative references (`Patient/mtb-<mrn>`).
- **Tag** every resource `http://elimu.io/fhir/tags#mtb-demo` (`mcode_helpers.add` does this).
- **A `display` on every coding.** The app shows bare codes otherwise. Use the official display (the validator rejects others).
- **One coherent timeline.** Dates are shifted into the demo period consistently (`mcode_helpers.shift_dates`); every date in the bundle, including those inside reports' text, must agree.
- **Validate** against `hl7.fhir.us.mcode#4.0.0` (0 errors expected; known exception: tx.fhir.org rejects frameshift p.HGVS), then load to `https://hapi.fhir.org/baseR4` and record the change in the scenario file.
- **Build with `mcode_helpers.py`** (`add`, `cc`, `qty`, `ext`, `write_bundle`, `relink`, `shift_dates`).

### 2.1 Genetic reports

Each genetic report is a structured report, not text only:

| Resource | Profile | Content |
|---|---|---|
| Genomics report | `mcode-genomics-report` (DiagnosticReport, LOINC 51969-4, categories Pathology and v2-0074 `GE`) | `result` references the Observations below; `specimen` references the tested specimen; `conclusion` is the lab's summary and must match the structured results exactly |
| Region studied | `mcode-genomic-region-studied` (LOINC 53041-0) | The genes tested (component 48018-6, `HGNC:` codes). Lets the app tell "tested and negative" from "not tested" |
| Genomic variant | `mcode-genomic-variant` (LOINC 69548-6, value LA9633-4 Present) | Gene (`HGNC:` form), c.HGVS with MANE transcript, p.HGVS with protein accession, SPDI in 81252-9 (same build as the server data), reference sequence, alleles, coordinate system, start, genomic source class (germline LA6683-2 / somatic LA6684-0), and VAF where the test reports one |
| Therapeutic implication | Genomics Reporting `therapeutic-implication` | In the server's CIViC shape: phenotype (Disease Ontology), level of evidence (CIViC text, e.g. "B (Clinical evidence)"), drug (RxNorm, capitalized name), predicted implication (CIViC text); `derivedFrom` the variant |

- **Which variants a report includes:** the reportable genes (the 228 genes in [Appendix: reportable genes](#appendix-reportable-genes)), as a real panel would report them.
- **Every implication must have a source.** Use a real CIViC (or label) evidence item and record it in the scenario file. Don't make one up. (HCC1143's FGFR2 pemigatinib implication was made up for the demo.)
- Only drug implications are carried; trial and diagnostic implications are left out.

### 2.2 Clinical content worth including

These made the syntheses harder or vaguer when missing:

- **A referral or visit reason** for the tumor board (e.g. a ServiceRequest with `reasonCode` or a note, or an Encounter `reasonCode`), so the reason for presenting the case is in the record, not inferred.
- **Values that decide treatment:** menopausal status (HR-positive breast cancer), Ki-67, PD-L1 on the most recent biopsy, current labs (liver, kidney, counts), ECOG, allergies with reactions, and enough treatment detail (cycles, dates, end dates).
- **Consistent staging:** the T category must fit the recorded size (MRN1234 records pT3 for a 2.5 cm tumor).
- **Germline testing** recorded whenever a variant could be inherited (L2345's BRCA2 Q66\* is a known germline pathogenic variant, reported as somatic, with no germline test).

## 3. The genomic data on the server

Loaded with the scripts in `utilities/` (see `utilities/README.md`):

1. **Molecular consequences:** run each VCF through `vcfPrepper` (the snpEff pipeline) first. Variants added by hand skip this and have no molecular consequences on the server.
2. **Variants:** `run_vcf2json` with one row per VCF in `vcfData.csv`: `vcf_filename`, `ref_build` (GRCh37 or GRCh38), `patient_id` (= MRN), `test_date` (yyyy-mm-dd), `test_id`, `specimen_id`, `genomic_source_class` (germline, somatic or mixed), `ratio_ad_dp`, `sample_position`. Outputs `variantsData.json`, `phaseData.json`, `molecularConsequences.json`.
3. **Tests:** one entry per test (`patientID`, `testDate`, `testID`, `specimenID`, `genomicBuild`, `studiedRegion` (a BED file, loaded with `bed2json`) and `dnaChangeType`). `run_vcf2json` doesn't update Patients or Tests. Watch for stray spaces in `dnaChangeType` values (the example file has `" INV"`, `" CNV"`).
4. Load to MongoDB (`to_mongo.py`).

`$find-study-metadata` returns these test fields, and with `ranges`, the studied and uncallable regions. The syntheses rely on them to say whether a negative finding means "absent" or "not tested".

## 4. Consistency between the report and the genomic data

These are the problems found in the four existing cases. A realistic case needs all of them right.

1. **Same test.** Each genetic report's test, specimen and date match a test on the server (test ID, `specimenID`, `testDate`). Found: HCC1143's server tests are dated 2020 against a 2024 report; L2345's single 2023 test holds the 2026 resistance mutations; MRN1234's test is dated 2023, before her 2026 diagnosis, from a specimen ("A2") not in the record. If a patient has two profiling reports (e.g. diagnosis and progression), load two tests.
2. **Variant types.** The test's `dnaChangeType` lists exactly the variant types in its data, and the report includes every reportable variant of those types, or states its reporting thresholds (e.g. "copy number ≥ 6 reported"). Found: L2345's test lists CNV and holds MET at 5 copies, but the report (curated from simple variants only) leaves it out. If the test only covers simple variants, the server shouldn't hold copy-number or structural calls for it.
3. **Regions.** `studiedRegion` matches what the test covers: a panel's genes for a panel, the genome for WGS. Found: MRN1234's report is a seven-gene germline panel, but its server test holds 17 somatic variants in other genes.
4. **Germline vs. somatic.** Germline and tumor results come from separate tests on separate specimens (blood vs. tumor), each with the right `genomic_source_class`. Every variant record carries its origin. Found: L2345's server variants have none.
5. **One genome build per patient**, where possible. Mixed builds force the server to lift ranges over between builds, and some ranges can't be (HCC1143: NCOA4, part of CRLF2 and PDE4DIP unsearchable). If builds must differ, record why.
6. **Realistic values:**
   - VAFs that vary between variants (not a fixed 0.25 for all, and not missing)
   - copy number on every copy-number and structural call
   - no Y-chromosome calls for a female patient (HCC1143 has 11)
   - no artifacts copied between patients (HCC1143 and L2345 share an identical TP53 3′ UTR deletion)
   - no duplicate representations of one variant (L2345 has L858R as both `T:G` and `TG:GT`)
   - SPDI coordinates that match the stated reference sequence (HCC1143 has a GRCh37 accession with GRCh38 coordinates)
7. **Report text and structure agree.** The `conclusion` names exactly the variants and implications in `result`, and the scenario file explains any deliberate difference between the report and the server data (e.g. a serious variant left off the report on purpose, so only the Genomics Explorer shows it).

## 5. After building a case

- Validate and load the bundle; load the genomic data; check the patient appears in the app's MRN box and that the Case Summary and Genomics Explorer agree with each other (or differ only as the scenario file says).
- Record decisions and deliberate quirks in `<MRN>-scenario.md`. The scenario file is for builders only: syntheses don't use it, because it wouldn't exist for a real patient.
- Write the case synthesis following `../public/data/synthesis/README.md`.

## Appendix: reportable genes

The 228 genes whose variants a genetic report in these cases includes (as given by Bob, 2026-11-05). A clinically serious variant in another gene is left off the report on purpose, so it can be found only in the Genomics Explorer.

ABL1, ACVR1, AIP, AKT1, AKT2, AKT3, ALK, AMER1, APC, AR, ARAF, ARID1A, ARID2, ASXL1, ATM, ATRX, AXIN2, B2M, BAP1, BARD1, BCL2, BCL9, BCOR, BLM, BMPR1A, BRAF, BRCA1, BRCA2, BRIP1, BTK, CALR, CARD11, CBFB, CCND1, CCND2, CCND3, CD79B, CDC73, CDH1, CDK12, CDK4, CDKN1B, CDKN2A, CHEK1, CHEK2, CIC, CREBBP, CSF1R, CTCF, CTNNA1, CTNNB1, CXCR4, CYLD, DDR2, DICER1, DNMT3A, EGFR, EP300, EPHA2, ERBB2, ERBB3, ERBB4, ERCC2, ESR1, EZH2, FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, FAS, FAT1, FBXW7, FGFR1, FGFR2, FGFR3, FGFR4, FH, FLCN, FLT1, FLT3, FOXA1, FOXL2, FUBP1, GATA3, GNA11, GNA13, GNAQ, GNAS, H3-3A, H3-3B, H3C2, HNF1A, HOXB13, HRAS, IDH1, IDH2, IRF4, JAK1, JAK2, JAK3, KDM5C, KDM6A, KDR, KEAP1, KIT, KLF4, KMT2A, KMT2C, KMT2D, KRAS, LZTR1, MAP2K1, MAP2K2, MAP2K4, MAP3K1, MAPK1, MAX, MED12, MEF2B, MEN1, MET, MITF, MLH1, MLH3, MPL, MRE11, MSH2, MSH3, MSH6, MTOR, MUTYH, MYC, MYCN, MYD88, NBN, NF1, NF2, NFE2L2, NFKBIA, NOTCH1, NPM1, NRAS, NSD1, NSD2, NTHL1, NTRK1, NTRK2, NTRK3, PALB2, PBRM1, PDGFRA, PDGFRB, PIK3CA, PIK3CB, PIK3R1, PIK3R2, PIM1, PMS2, POLD1, POLE, POT1, PPP2R1A, PPP2R2A, PRDM1, PRKACA, PRKAR1A, PRKDC, PTCH1, PTEN, PTPN11, RAC1, RAD50, RAD51B, RAD51C, RAD51D, RAD54L, RAF1, RASA1, RB1, RET, RHOA, RNF43, ROS1, RUNX1, SDHA, SDHAF2, SDHB, SDHC, SDHD, SETD2, SF3B1, SMAD2, SMAD4, SMARCA4, SMARCB1, SMARCE1, SMO, SOCS1, SOS1, SPEN, SPOP, SRC, STAG2, STAT3, STK11, SUFU, TCF7L2, TERT, TET2, TMEM127, TNFAIP3, TNFRSF14, TP53, TRAF7, TSC1, TSC2, U2AF1, VHL, WRN, WT1, XPO1, XRCC1
