# MTB app terminology

How the MTB app uses clinical terminologies to:
- decide whether a genomic implication is actionable for the patient's tumor type
- choose gene-list presets
- tell the oncogenicity predictor the tumor type

## Overview

```
Patient's primary cancer Condition (SNOMED CT code)          or a generic tumor type chosen in the Dx dropdown
  │
  │  MTB_TumorTypeCrosswalk.csv  (built from the VSAC value sets in terminology/value-sets/)
  ▼
Tumor type: Breast Cancer | Non-Small Cell Lung Cancer | Hepatocellular Carcinoma | Colorectal Cancer
  │
  ├─ Presets ─────────────── MTB_KB_GeneLists.csv rows with that Cancer category
  ├─ Oncogenicity predictor ─ the crosswalk's predictorTumorType, or nothing
  └─ Actionability ───────── MTB_TumorTypeImplicationCodes.csv: the Disease Ontology codes
                              allowed for that tumor type, plus the pan-cancer codes
```

The **tumor type** is the key that ties everything together. It is the `geneListCategory` in the crosswalk and the implication codes, and the `Cancer category` in the gene lists. All three files must spell it identically.

## 1. Tumor type and value sets

### From diagnosis to tumor type

The app looks up the patient's primary cancer Condition in `MTB_TumorTypeCrosswalk.csv` (`mapCodingsToTumorType`, `src/cancerFilter/tumorTypeCrosswalk.ts`). The first coding whose system and code **exactly** match a row gives the tumor type.

- A code under a listed code doesn't count, so a diagnosis coded more specifically than the value set gets no tumor type.
- Only SNOMED CT is mapped. ICD-10-CM and ICD-O-3 codes on the Condition are ignored here (see [Future topics](#future-topics)).
- The user can override the tumor type with a generic entry in the Dx dropdown. There is one per tumor type in the crosswalk.
- With no tumor type, the presets offer only the non-cancer lists, the predictor gets no tumor type, and "Actionable, this tumor type" isn't available (see section 3).

### Files in `terminology/value-sets/`

Each file is a VSAC expansion spreadsheet, used exactly as downloaded: a "Value Set Info" sheet with the OID, and an "Expansion List" sheet with the codes. `build_crosswalk.py` lists each file with the tumor type and predictor tumor type it maps to.

| File | VSAC value set (steward) | OID | Codes | Tumor type | Predictor tumor type |
|---|---|---|---|---|---|
| `BreastCancer.xlsx` | Breast Cancer | 2.16.840.1.113883.3.1444.3.102 | 29 SNOMED CT | Breast Cancer | Breast Cancer |
| `NonSmallCellLungCancer.xlsx` | Non-Small Cell Lung Cancer | 2.16.840.1.113883.3.1444.3.277 | 89 SNOMED CT | Non-Small Cell Lung Cancer | Non-Small Cell Lung Cancer |
| `HepatocellularCarcinoma.xlsx` | Hepatocellular Carcinoma | 2.16.840.1.113762.1.4.1056.93 | 5 SNOMED CT | Hepatocellular Carcinoma | *(none)* |
| `ColorectalCancer.xlsx` | Colorectal Cancer (ASCO) | 2.16.840.1.113883.3.1444.3.215 | 22 SNOMED CT | Colorectal Cancer | *(none)* |

**Local codes.** The colorectal value set has general neoplasm codes but no adenocarcinoma codes. VSAC's broader colorectal set also takes in cancers that don't apply, such as colon lymphoma. So `build_crosswalk.py` adds three SNOMED CT codes of its own to Colorectal Cancer (its `LOCAL_CODES` list), with `valueSetOid` set to `local`:

| SNOMED CT | Display |
|---|---|
| 1701000119104 | Primary adenocarcinoma of colon |
| 254582000 | Adenocarcinoma of rectum |
| 425178004 | Adenocarcinoma of rectosigmoid junction |

**Rebuilding.** The build fails if a code appears in two value sets, or in a value set and the local list. When that happens because a newer VSAC version added a local code, remove the code from the local list. To refresh:
1. Download newer expansions from VSAC.
2. Replace the spreadsheets.
3. Run:

```
python terminology/build_crosswalk.py
```

## 2. Spreadsheets in `public/data/`

| File | Built by | What it's for | Columns |
|---|---|---|---|
| `MTB_TumorTypeCrosswalk.csv` | `terminology/build_crosswalk.py` (don't edit by hand) | Diagnosis code → tumor type and predictor tumor type | `system, code, display, predictorTumorType, geneListCategory, valueSetOid` |
| `MTB_TumorTypeImplicationCodes.csv` | By hand, from CIViC's disease list | The cancer-type codes that make an implication match each tumor type (section 3) | `geneListCategory, system, code, display` |
| `MTB_KB_GeneLists.csv` | By hand | Preset gene lists per tumor type (section 4) | `Cancer category, Label, Region` |
| `CIViC_Cooccurrence_KB.csv` | `scripts/build_civic_cooccurrence_kb.py` | CIViC evidence that needs two or more variants together | Includes `disease` as a **name only**, with no code. Co-occurrence evidence isn't filtered by tumor type. |

## 3. Actionability: "Actionable, this tumor type"

### What an implication carries

Every therapeutic or diagnostic implication codes the cancer types its evidence applies to. There is one LOINC `81259-4` component per cancer type, and an implication can have several.

| Source | Coding |
|---|---|
| Operations (`$find-subject-tx-implications`, `$find-subject-dx-implications`), from CIViC | `https://disease-ontology.org`, e.g. `DOID:3908` Lung Non-small Cell Carcinoma |
| The patient's genetic reports | The same Disease Ontology codes. Solid Tumor is `http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl` `C9292` Solid Neoplasm. |

### The match

An implication is actionable for this tumor type when both of these hold:

1. **Its evidence level is A:** the text starts with "A ", equals "A", or contains "Validated association" (`isALevelEvidence`).
2. **Any of its cancer-type codings** (system and code) is in `MTB_TumorTypeImplicationCodes.csv`, either for the selected tumor type or for `Any` (`matchesTumorType`, `src/cancerFilter/tumorTypeCodes.ts`).

A variant is shown when at least one of its implications passes. An implication with no coded cancer type never passes this filter. It still shows under "Actionable, any tumor type".

**Leading zeros.** Disease Ontology codes are compared without leading zeros in the number. The genomics server sends e.g. `DOID:50861` where CIViC, the reports and the CSV have `DOID:0050861`.

**No tumor type.** With no tumor type, "Actionable, this tumor type" is greyed out in the Filters sidebar and "Actionable, any tumor type" applies instead. The choice is kept, and applies again once the Dx has a tumor type.

**Where it applies.** The filter acts on the Genomics Explorer only. The Case Summary's reported variants card shows every implication, unfiltered.

**Pan-cancer codes.** These are allowed for every tumor type (category `Any`). An example is an A-level immunotherapy approved for all solid tumors.

| System | Code | Display | CIViC disease |
|---|---|---|---|
| `https://disease-ontology.org` | DOID:162 | Cancer | Cancer |
| `http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl` | C9292 | Solid Neoplasm | Solid Tumor (CIViC gives it no DOID) |

**Per-tumor-type codes.** The codes are taken from CIViC's disease list.

| Tumor type | Disease Ontology codes |
|---|---|
| Colorectal Cancer | DOID:9256 Colorectal Cancer, 0050861 Colorectal Adenocarcinoma, 219 Colon Cancer, 234 Colon Adenocarcinoma, 1520 Colon Carcinoma, 3029 Colon Mucinous Adenocarcinoma, 3033 Colon Signet Ring Adenocarcinoma, 1993 Rectum Cancer, 1996 Rectum Adenocarcinoma |
| Non-Small Cell Lung Cancer | DOID:3908 Lung Non-small Cell Carcinoma, 3910 Lung Adenocarcinoma, 3907 Lung Squamous Cell Carcinoma, 6482 Lung Acinar Adenocarcinoma, 0080303 Mucinous Lung Adenocarcinoma, 4556 Lung Large Cell Carcinoma, 4829 Adenosquamous Lung Carcinoma, 0080777 Lung Sarcomatoid Carcinoma, 0050932 Lung Mucoepidermoid Carcinoma; broader: 1324 Lung Cancer, 3905 Lung Carcinoma |
| Breast Cancer | DOID:1612 Breast Cancer, 3459 Breast Carcinoma, 3458 Breast Adenocarcinoma, 3007 Breast Ductal Carcinoma, 0050938 Breast Lobular Carcinoma, 4922 Breast Secretory Carcinoma, 6587 Breast Tubular Carcinoma, 6263 Inflammatory Breast Carcinoma, 0060548 Luminal Breast Carcinoma A, 0060075 ER-positive, 0060076 ER-negative, 0060079 HER2-positive, 0060080 HER2-negative, 0060081 Triple-negative |
| Hepatocellular Carcinoma | DOID:684 Hepatocellular Carcinoma, 0080182 Mixed Fibrolamellar Hepatocellular Carcinoma, 0070322 Childhood Hepatocellular Carcinoma; broader: 3571 Liver Cancer |

Choices made in building these lists:
- **Broader terms are included for NSCLC and HCC.** CIViC evidence filed under "Lung Cancer" or "Liver Cancer" is usually imprecise wording for the same disease.
- **Excluded:** small cell lung types (DOID:5409, 0050917), Breast Angiosarcoma (DOID:4511, a sarcoma treated very differently), Hepatoblastoma (DOID:687) and Intrahepatic Cholangiocarcinoma (DOID:4928).
- **Breast receptor status isn't checked.** Receptor-status subtypes match any breast patient.

Keep this table in step with the CSV.

**Matching is exact.** It doesn't walk the Disease Ontology hierarchy. A CIViC disease that should count for a tumor type must be listed explicitly. When CIViC adds a disease, check whether it belongs in a list.

## 4. Presets

`MTB_KB_GeneLists.csv` (`src/cancerFilter/geneListKb.ts`) has one row per gene. The preset buttons offer the lists whose `Cancer category` equals the tumor type, plus the `Any` lists, which are offered whatever the tumor type.

| Cancer category | Lists (genes) |
|---|---|
| Breast Cancer | Actionable Gene List (10), Extended Gene List (17) |
| Non-Small Cell Lung Cancer | Actionable Gene List (11), Extended Gene List (16) |
| Hepatocellular Carcinoma | Extended Gene List (25) |
| Colorectal Cancer | Actionable Gene List (8), Extended Gene List (108) |
| Any | ACMG Secondary Findings (83) |

The colorectal lists were built from CIViC evidence (accepted and submitted) on the colorectal Disease Ontology codes above:
- **Actionable:** genes with level-A evidence (BRAF, ERBB2, KRAS, NRAS), plus the mismatch-repair genes MLH1, MSH2, MSH6 and PMS2. These stand in for MSI-High, which has level-A evidence but is a biomarker, not a gene.
- **Extended:** genes with any evidence, plus everything in Actionable.
- **Fusions:** only the kinase partner is listed.

## 5. Oncogenicity predictor

The crosswalk's `predictorTumorType` is sent as the predictor's optional `tumorType` (`src/lib/oncogenicity.ts`). When it's blank, because the predictor has no matching type, **no tumor type is sent**, and the predictor uses its tumor-type-independent rules. Predictions are cleared when the Dx changes.

The predictor accepts these values, from its OpenAPI description:

> Bladder Carcinoma, Breast Cancer, Granulosa Cell Tumor, Hairy Cell Leukemia, Hodgkin Lymphoma, Neuroblastoma, Non-Small Cell Lung Cancer, Parathyroid Carcinoma, Peripheral T-Cell Lymphoma, Retinoblastoma, Renal Cell Carcinoma, T-Cell Acute Lymphoblastic Leukemia, T-Cell Lymphoblastic Lymphoma, Urothelial Carcinoma

Of the app's tumor types, only Breast Cancer and Non-Small Cell Lung Cancer are sent. Hepatocellular Carcinoma and Colorectal Cancer send none. If the predictor adds a type, set it for that value set in `build_crosswalk.py` and rebuild.

## Future topics

### Knowing whether a code is a cancer code

A Condition is treated as a primary cancer today when any of these is true (`getPrimaryCancerConditions`, `src/lib/clinicalSummary.ts`):
- It has the mCODE `mcode-primary-cancer-condition` profile.
- Its code is in the crosswalk.
- It has an ICD-10-CM code in C00–C76 or C80–C96. This is a pattern, not a value set.

Secondary cancers are recognized by the mCODE `mcode-secondary-cancer-condition` profile or ICD-10-CM C77–C79.

**Gap:** a SNOMED-only primary cancer outside our four tumor types, without the mCODE profile, isn't recognized. A general "primary malignant neoplasm" value set, for SNOMED CT and ICD-10-CM, would close it. Recognizing a cancer is separate from mapping it to a tumor type: a recognized cancer can still have no tumor type.

### Other mCODE code systems: ICD-10-CM and ICD-O-3

mCODE codes a primary cancer as follows:
- `Condition.code`: SNOMED CT or ICD-10-CM
- Histology (the `mcode-histology-morphology-behavior` extension): SNOMED CT or ICD-O-3 morphology
- `bodySite`: SNOMED CT or ICD-O-3 topography

Our patients carry these alongside SNOMED CT, for example:

| Patient | Condition.code | Histology |
|---|---|---|
| TCGA colorectal (e.g. TCGA-A6-2674) | SNOMED 1701000119104; ICD-10-CM C18.7 | ICD-O-3 8480/3 Mucinous adenocarcinoma |
| TCGA-DD-A1EH | SNOMED 109841003; ICD-10-CM C22.0 | SNOMED 25370001; ICD-O-3 8170/3 Hepatocellular carcinoma, NOS |
| HCC1143 | SNOMED 408643008; ICD-10-CM C50.411 | SNOMED 82711006; ICD-O-3 8500/3 |

Ways to extend:
- **ICD-10-CM in the crosswalk.** Many VSAC value sets have ICD-10-CM versions, or include both code systems. Add the system to `CODE_SYSTEMS` in `build_crosswalk.py` and download those expansions. ICD-10-CM can't tell histology apart, though: C34.x covers both small cell and non-small cell lung cancer. So it fits tumor types defined by site (Breast, Colorectal) better than those defined by histology (NSCLC, HCC).
- **ICD-O-3 topography plus morphology.** Site and histology together pin down a tumor type, for example C22.0 with 8170/3 for HCC. That needs a two-part rule rather than a single-code value set.
- **Order of precedence.** Decide which wins when a Condition's codes map to different tumor types. Today, the first matching SNOMED CT coding wins.

The Operations server can translate SNOMED CT and ICD-10 conditions to Disease Ontology codes (its `conditions` parameter, via concept maps). The app doesn't use this. It could someday give a direct diagnosis-to-DOID path, but exact DOID matching would still need the hierarchy, or curated lists, to connect a specific diagnosis to broader evidence.
