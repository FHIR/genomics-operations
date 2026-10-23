# Jenny M (MRN1234): clinical scenario

Synthetic patient for the MTB app clinical-data demo (see `mtb-clinical-data-design.md` in the repo root). Jenny M is the mCODE IG's own illustrative patient.

**Files**
- [MRN1234-transaction-bundle.json](MRN1234-transaction-bundle.json): the bundle (mCODE STU4 / US Core 6.1), **edited directly**.
  It was produced once from Bob's draft (assembled from the IG's Jenny M examples) by a transform script applying the revisions below. The draft and script were then retired (2026-10-22), so this JSON is now the source of truth.

## Base scenario (mCODE IG)

From the [mCODE IG examples page](https://build.fhir.org/ig/HL7/fhir-mCODE-ig/examples.html):

> Patient Jenny M. is a 55 year old non-Hispanic white female.

- **History:** depression, anxiety, hypertension; current smoker, 20 pack-years. Maternal aunt with ovarian cancer (69), sister with breast cancer (64), paternal uncle died of pancreatic cancer.
- **Diagnosis:** screening mammogram in February 2018 showed a possible mass. Ultrasound-guided biopsy found invasive ductal carcinoma, ER+, PR−, HER2−. Clinical cT3 cN0 cM0. A seven-gene panel found a pathogenic PALB2 variant (c.3549C>A). Oncotype DX recurrence score 47.
- **Surgery:** partial mastectomy: 2.5 cm, grade 2, no lymphovascular invasion, margins negative, 0 of 3 sentinel nodes. Pathologic pT3 pN0, stage IIB. ECOG 0.
- **Treatment:** dose-dense AC (doxorubicin 60 mg/m², cyclophosphamide 600 mg/m²) × 4, then paclitaxel 175 mg/m². Whole-breast radiation with regional nodal irradiation.
- **Follow-up:** one month after completing first-line therapy, disease status improving; plan surveillance imaging.

## Timeline (dates moved to the present)

As with TCGA-DD-A1EH, every date moves by the same 2,905 days (birth date included), so intervals and age are unchanged. The last data point, 2018-11-01, becomes 2026-10-15.

| Bundle date | Draft/IG date | Event |
|---|---|---|
| 1972-12-15 | 1965-01-01 | Date of birth |
| 2026-01-15 | 2018-02-01 | Screening mammogram |
| 2026-02-17 | 2018-03-06 | Ultrasound-guided breast biopsy. Weight |
| 2026-02-21 | 2018-03-10 | ER positive, PR negative, HER2 negative |
| 2026-02-26 | 2018-03-15 | Seven-gene panel: PALB2 c.3549C>A |
| 2026-02-27 | 2018-03-16 | Diagnosis; clinical cT3 cN0 cM0, IIB. Smoking status, comorbidities |
| 2026-03-15 | 2018-04-01 | Partial mastectomy + sentinel nodes: 2.5 cm, grade 2, margins negative, 0/3 nodes |
| 2026-03-19 | 2018-04-05 | Pathologic stage pT3 pN0, IIB |
| 2026-03-26 | 2018-04-12 | ECOG 0; Oncotype DX 47. Dose-dense AC starts |
| 2026-05-03 | 2018-05-20 | Neutrophils 3,000/µL |
| 2026-05-07 | 2018-05-24 | Last AC cycle |
| 2026-05-21 – 07-02 | 2018-06-07 – 07-19 | Paclitaxel, dose-dense × 4 |
| 2026-07-29 – 10-08 | 2018-08-15 – 10-25 | Radiotherapy: left whole breast 60 Gy / 30 fractions, regional nodes 50 Gy / 25 fractions |
| 2026-09-13 | 2018-09-30 | Anastrozole starts (active) |
| 2026-10-15 | 2018-11-01 | Disease status: improved |

## Revisions made by the build script

The draft failed validation with 69 errors. 53 of them came from one cause: the Patient's identifier used an example URL, so it failed to match US Core and every reference to it failed too. The draft also used POST with server ids and absolute references to the authoring server, and it referenced 15 resources at `http://example.org/fhir/...` that weren't in the bundle.

**Structure (matching the other demo bundles)**
- Fixed `PUT` ids prefixed `mtb-mrn1234`, relative references, `mtb-demo` tag, server metadata removed. MRN system changed from the example URL `http://hospital.example.org` to `http://elimu.io/fhir/mrn`.
- The duplicate Dr. Owen Oncologist (two identical Practitioners) is merged into one.

**Resources the draft referenced but didn't contain**, added with values from the IG's Jenny M example:
- Practitioners Mary Obgyn, Jane Radiotech and Kyle Anydoc; Organization Bedrock Medicine.
- The surgical specimen; the Nottingham grade 2 observation.
- Depression, anxiety and hypertension Conditions (the comorbidities observation pointed at them).
- The pathologic stage: pT3, pN0 and stage group IIB. The IG's pM was a data-absent-reason code inside `valueCodeableConcept`, so the pathologic group uses the clinical cM0 instead.
- The two radiotherapy volumes, now **left whole breast** and **left regional nodes**. The IG had chest wall and mediastinal nodes, which don't fit a partial mastectomy.
- The sister's breast cancer (age 64), which is in the narrative but was missing from the draft.

**Content fixes**
- **Genomic resources dropped:** GenomicVariant (PALB2), genomic region studied and genomics report, replaced by a text-only genomic report (design doc §3.3). See the genomics-server note below.
- **Primary cancer:**
  - The SNOMED code from the IG, 353431000119107, is now inactive. It's replaced by 408643008 *Infiltrating duct carcinoma of breast*, a code in the VSAC Breast Cancer value set that sets the predictor's tumor type (design doc §5.5). Laterality is in the text ("…of left breast") and the body site.
  - ICD-10-CM C50.912 and the left-breast body site are added.
  - The stage summary said **IIIC**, contradicting the stage group (IIB). Now there is one IIB entry per stage group.
- **Clinical stage group:** its members pointed at the IG's *pathologic* T/N/M; now the clinical cT3/cN0/cM0. Category, focus and the current AJCC method display are added.
- **Chemotherapy:** the draft had all three drugs `active` and starting together (doxorubicin and cyclophosphamide on the surgery date). Now dose-dense AC × 4 (completed), then paclitaxel × 4 (completed), per the narrative. The doxorubicin dose moves from a rate to a dose.
- **Radiotherapy:** the reason was ICD-10 C50.811 (*right* breast) and the site chest wall. Now the reason is the left-breast primary and the site is the left breast.
- **Neutrophils:** 3000 with unit 10³/µL is 1,000× too high; the unit is now /µL.
- **Other:** ECOG and comorbidities categories; current displays (ovarian and pancreatic cancer, HER2); tumor markers linked to the cancer and the biopsy specimen; a pathology report conclusion.

## Inconsistencies kept from the IG (flagged, not changed)

- **Size vs. T category:** the tumor is 2.5 cm, which is T2 (2–5 cm), yet the IG says cT3 and pT3. With pT2 pN0 the pathologic stage would be IIA. The bundle keeps the IG's T3 / IIB.
- **Age:** the narrative says 55, but the IG's birth date (1965-01-01) makes her 53 at diagnosis. The bundle keeps the birth date (shifted), so she's 53.
- **PALB2 origin:** the IG marks the variant *somatic*. Seven-gene panels like this are usually germline tests. The text report doesn't state an origin.

## Genomic data on the genomics operations server

Queried 2026-10-22 (subject `MRN1234`). The server has data for this subject, but **not the PALB2 variant**: the only hit across PALB2, BRCA1 and TP53 is TP53 NC_000017.10:7578412:C:A (GRCh37). So the server's MRN1234 is a different synthetic profile from the IG's Jenny M. The MTB demo would show a TP53 variant in Genomic findings and a PALB2 finding in the clinical reports. That's worth reconciling, or at least explaining, before a demo.

## Data-shape features this bundle exercises

- **Pathology detail Observations** (invasion, margins, sentinel nodes, DCIS) under a pathology DiagnosticReport.
- **Oncotype DX** coded with a GTR code (not LOINC). It will fall into Labs by category.
- **Three FamilyMemberHistory** resources, an Organization performer, and a pack-years observation.
- **Primary cancer `clinicalStatus` = remission** with a CancerDiseaseStatus of improved.
- **An active endocrine therapy** (anastrozole) after completed chemotherapy and radiation.

## Corrections (2026-10-23)

- **Missed dates:** the original date shift skipped values that include a time of day (e.g. `…T00:00:00Z`). 10 such values were still in the original years; now shifted by the same 2,905 days.
- **Missing displays:** 27 codings had a code but no `display` (e.g. LOINC 89247-1 ECOG, LOINC 33728-7 tumor size dimension), so the app fell back to showing the bare code. Displays added from tx.fhir.org `$lookup`. CPT 99213 on the Encounters was left without one (licensed text; the Encounters carry `text`).
- **Sentinel nodes:** the observation *Sentinel lymph nodes with metastasis [#]* had no value, only a component for nodes examined (3). Added its value: 0 (the scenario's 0 of 3), keeping the examined count as the component.
- Revalidated (0 errors) and reloaded to hapi.fhir.org.

## Validation

HL7 FHIR validator against `hl7.fhir.us.mcode#4.0.0` with tx.fhir.org: **0 errors** (2026-10-22; the draft had 69). The remaining warnings are the benign kinds seen in the other bundles: no narrative, RxNorm version mismatch (the draft uses RxNorm ingredient codes), extensible value-set notes, and the GTR code system being unknown to the validator.

**Confirmed approach:** follows Bob's decisions for TCGA-DD-A1EH (dates moved to the present, genomic resources off HAPI, shared MRN system).

## Loading (loaded to hapi.fhir.org 2026-10-22)

```bash
curl -X POST -H "Content-Type: application/fhir+json" \
  --data @MRN1234-transaction-bundle.json https://hapi.fhir.org/baseR4
```
