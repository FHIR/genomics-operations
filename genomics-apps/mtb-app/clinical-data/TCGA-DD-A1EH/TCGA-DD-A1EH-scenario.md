# Higado Sobreviviente (MRN TCGA-DD-A1EH): clinical scenario

Synthetic patient for the MTB app clinical-data demo (see `mtb-clinical-data-design.md` in the repo root).

**Files**
- [TCGA-DD-A1EH-transaction-bundle.json](TCGA-DD-A1EH-transaction-bundle.json): the bundle (mCODE STU4 / US Core 6.1), **edited directly**.
  It was produced once from Bob's hand-authored bundle (2024) by a transform script applying the revisions below. The original and script were then retired (2026-10-22), so this JSON is now the source of truth.

## Base scenario (from Bob's case description)

> This case is derived from [TCGA-DD-A1EH](https://www.cbioportal.org/patient/summary?studyId=lihc_tcga_pan_can_atlas_2018&caseId=TCGA-DD-A1EH), modified for illustrative purposes. Patient Higado Sobreviviente (DOB 1998-04-14) is a Hispanic male with a known history of chronic Hep B, who was found in June 2021 at the age of 23 to have a liver mass, with biopsy showing hepatocellular carcinoma.

**Dates are moved to the present (Bob, 2026-10-22).** Every date shifts by the same 1,644 days, including the birth date, so intervals and the patient's age (23 at diagnosis) are unchanged. The last visit moves from 2022-04-15 to 2026-10-15. The case description keeps its original dates; the bundle uses the shifted ones:

| Bundle date | Case description date | Event |
|---|---|---|
| 2002-10-14 | 1998-04-14 | Date of birth |
| 2025-12-20 | 2021-06-20 | Initial presentation; liver mass noted and biopsied. Chronic hepatitis B; latex allergy |
| 2026-01-06 | 2021-07-07 | Initial diagnosis: HCC, clinical T3 N0 M0, stage IIIA. Child-Pugh 5 (class A), AFP, ECOG 0. Begin tenofovir |
| 2026-01-26 | 2021-07-27 | Surgical right lobe resection. Pathologic T3 N0 M0; inactive chronic hepatitis with fibrosis |
| 2026-05-02 | 2021-10-31 | Liver and lung recurrence. Begin sorafenib. Disease status: worsened |
| 2026-07-14 | 2022-01-12 | Clinical progressive disease. NGS ordered. Disease status: worsened |
| 2026-07-28 | 2022-01-26 | NGS results reviewed (MSS, TMB 1.3, fusion, simple variants); referred for clinical trial. ECOG 0 |
| 2026-10-15 | 2022-04-15 | Last visit: on trial NCT05023655, tazemetostat started, sorafenib stopped. Disease status: stable |

Hepatitis B onset (year only) moves from 2017 to 2022.

## Revisions made by the build script

The source bundle failed validation with 43 errors (20 in the genomic resources). It also used POST with server-assigned ids and absolute references to the authoring server, so loading it into HAPI would create new ids and duplicates on every reload. The script applies these changes:

**Structure (matching the other demo bundles)**
- Fixed `PUT` ids prefixed `mtb-tcga-dd-a1eh`, relative references, `meta` replaced (profile kept, `mtb-demo` tag added, server versionId/lastUpdated/source removed).
- MRN identifier system changed from `http://elimu.io/gacs` to the shared `http://elimu.io/fhir/mrn`. Patient profile is mCODE CancerPatient (plus birth sex).

**Confirmed by Bob 2026-10-22:** dates moved to the present, genomic resources dropped from the HAPI bundle, the two primary cancers merged, the MRN system changed.

**Content**
- **Genomic resources dropped (Bob, 2026-10-22):** 36 GenomicVariant, 3 therapeutic-implication and 3 diagnostic-implication Observations. The genomics server already holds this patient's variants (ARID1A, AXIN1, SI, PACS2 and others, on GRCh37 accessions), and genomics isn't duplicated on HAPI for the demo. Real-world feeds *will* carry genomic resources alongside clinical data; how the app handles that overlap is an open design item (design doc §3.3). MSI and TMB are kept as tumor markers. A **text-only NGS DiagnosticReport** summarizes the findings and candidate trials. It lists only variants the genomics server holds (ARID1A p.Arg1721* and p.Val2169fs; AXIN1 p.Arg647fs).
  - Revised 2026-10-22: the source bundle's APOB::CLU fusion isn't on the genomics server, and its second "AXIN1" variant (c.715C>T) was placed on chromosome 17, not at AXIN1 (chromosome 16). Both were removed from the report.
- **One primary cancer:** the source had two PrimaryCancerConditions, the original (right lobe, stage IIIA) and a recurrence copy (left lobe, onset 2021-10-31, stage IV, category encounter-diagnosis). The app would read two separate primary cancers. Now there is one primary with `clinicalStatus` = recurrence. The left-lobe recurrence is a **SecondaryCancerCondition** (metastatic malignant neoplasm to liver, left lobe). The stage-IV summary is dropped: AJCC stage isn't reassigned at recurrence.
- **Primary cancer details added:** histology (SNOMED 25370001 Hepatocellular carcinoma; ICD-O-3 8170/3), ICD-10-CM C22.0, verification status, links to both stage groups.
- **TNM:** added the missing `category` (9 errors), `focus` and `performer`; updated the AJCC 8th edition method display (the old one is inactive). The pathologic M (a data-absent-reason code placed in `valueCodeableConcept`) is replaced by the clinical cM0 in the pathologic stage group.
- **ECOG category:** the source used a US Core 3.x code system that no longer exists; now `survey`.
- **Disease status:** added category, evidence type (imaging) and performer.
- **Tumor markers:** linked to the cancer; TMB gets its unit (mutations/Mb).
- **Latex allergy:** category `environment` (was medication); verification status added.
- **Procedures:** performer, reason; resection gets curative intent and body site.
- **ResearchStudy:** the source's identifier was a placeholder (NCT01234567). Its title is that of **NCT05023655** (*Phase II Study of Tazemetostat in Solid Tumors Harboring an ARID1A Mutation*), so the identifier now matches.
- **Medications added** (in the case description, missing from the source):
  - tenofovir 300 mg daily (2021-07-07, active; US Core, for hepatitis B)
  - sorafenib 400 mg twice daily (2021-10-31 → 2022-04-15, stopped)
  - tazemetostat 800 mg twice daily (2022-04-15, active, on trial)

**Dates:** all shifted to the present by a fixed 1,644 days (see the table above). The shift is the last step of the build script, so the source bundle keeps its original dates.

**Left as authored:** the six Encounters, the Practitioner (Dr. Owen Oncologist), Child-Pugh (as a lab observation), AFP, the comorbidities observation, the ResearchSubject (start 2022-02-07, now 2026-08-09), and the Patient's race, ethnicity and address.

## Genomic data on the genomics operations server

Queried 2026-10-22 (subject `TCGA-DD-A1EH`). Variants returned include ARID1A NC_000001.10:27105549:C:T and NC_000001.10:27106893::T, AXIN1 NC_000016.9:347071:GG:GGG, SI NC_000003.11:164750328:G:A and PACS2 NC_000014.8:105858083:G:A. These are the same variants as the source bundle's GenomicVariant resources. Coordinates are on GRCh37 accessions.

## Trials (ClinicalTrials.gov, checked 2026-10-22)

| Trial | Title | Status today |
|---|---|---|
| NCT05023655 | Phase II Study of Tazemetostat in Solid Tumors Harboring an ARID1A Mutation | Terminated |
| NCT06005974 | A Study of REC-4881 in Participants With Cancers Which Have an AXIN1 or APC Mutation | Terminated |

The bundle keeps ResearchStudy status `active`. With the dates moved to 2026, this is illustrative: both trials have since been terminated in reality.

## Data-shape features this bundle exercises

- **Encounters, ResearchStudy, ResearchSubject:** resource types the other bundles don't have. They should land in "Other resources", or a future section.
- **LOINC staging codes** (21905-5, 21899-0, …) for TNM, where the other bundles use SNOMED.
- **Primary-cancer recurrence plus an intrahepatic SecondaryCancerCondition.**
- **A background liver disease (hepatitis B) recorded as a comorbidity**, with an active antiviral.

## Validation

HL7 FHIR validator against `hl7.fhir.us.mcode#4.0.0` with tx.fhir.org: **0 errors** (2026-10-22; the source had 43). The remaining warnings are the benign kinds seen in the other bundles: no narrative, RxNorm version mismatch, extensible value-set notes, partial ICD-O-3 support on tx.fhir.org, and the Child-Pugh interpretation code (SNOMED) falling outside the base interpretation value set.

## Loading (loaded to hapi.fhir.org 2026-10-22)

```bash
curl -X POST -H "Content-Type: application/fhir+json" \
  --data @TCGA-DD-A1EH-transaction-bundle.json https://hapi.fhir.org/baseR4
```
