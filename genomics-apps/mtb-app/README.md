# Molecular Tumor Board (MTB) — Variant & Implication Explorer

Live demo: [https://mtb-app-elimu1.vercel.app/](https://mtb-app-elimu1.vercel.app/)

A web application that helps clinicians and researchers prepare a patient for a Molecular Tumor Board (MTB): the clinical record, the lab's reported genomic findings, the patient's full genomic data with **therapeutic (Tx)**, **diagnostic (Dx)** and **molecular consequence** annotations, and an AI case synthesis. It's a proof of concept showing how MTB features derive from FHIR: clinical data from a FHIR server (mCODE), genomic data from the [Genomics Operations](https://github.com/FHIR/genomics-operations) reference implementation.

---

## Table of Contents

* [Overview](#overview)
* [Core Features](#core-features)
* [Data Sources](#data-sources)
* [Tech Stack](#tech-stack)
* [Project Structure](#project-structure)
* [Usage](#usage)
* [Development](#development)
* [License](#license)

---

## Overview

The app is built around a typical MTB workflow: a patient's case, including tumor genomic profiling, is presented, and the board reviews the findings, their therapeutic and diagnostic implications, and the treatment options.

* **Patient banner:** choose a patient by MRN. The banner shows name, age, sex, latest stage, ECOG, disease status and any study enrollment, and the Dx (set from the patient's recorded primary cancer, which can be overridden with a generic tumor type).
* **Four tabs:**
  * **Case Summary:** the clinical record from the FHIR server, including the **Reported variants** card: the variants and therapeutic implications from the lab's genetic reports.
  * **Genomics Explorer:** read-only search of the patient's complete genomic data set by gene symbol or genomic range, with presets, actionability filters, co-occurring variants, pathway diagrams and oncogenicity prediction. Reported variants are layered into the results.
  * **Imaging:** a placeholder for pathology slides and diagnostic images.
  * **AI Synthesis:** a case synthesis for the tumor board (pre-generated for the demo patients).

The Case Summary shows what the lab chose to report, as of the report date; the Genomics Explorer shows the complete data set with current annotations. A variant in the Explorer that isn't on the report is expected, not an error.

---

## Core Features

* **Patient Banner and Dx**: MRN selection (patients with clinical data listed first), key clinical status, and a Dx that drives the gene-list presets, the "Actionable, this tumor type" filter and the oncogenicity predictor's tumor type. Recorded diagnoses map to tumor types by SNOMED CT code (`public/data/MTB_TumorTypeCrosswalk.csv`)
* **Case Summary**: diagnosis and staging, treatment timeline, labs and pathology, reports and other clinical sections, with a jump bar
* **Reported Variants**: the lab's reported variants grouped by report, with origin, therapeutic implications and evidence levels, and a link to the full data in the Genomics Explorer
* **Reported Variant Layering**: reported variants in the searched genes or ranges join the Explorer's results (a variant in both appears once), marked with a Reported chip; report implications carry a Report tag and count toward the actionability filters. A Reported genes preset and a Reported variants actionability option
* **Intelligent Search**: Search by gene symbols (e.g., `BRAF`, `EGFR`) or genomic ranges (e.g., `NC_000007.14:55019016-55211628`)
* **Guided Search Presets**: Cancer-type-specific Actionable Gene List and Extended Gene List buttons populate the search box with predefined terms
* **Advanced Filtering**: Sidebar actionability, evidence levels, medications, implications, and molecular consequences
* **Real-time Results**: Incremental loading with search status indicators and cancellation support
* **Comprehensive Data Display**: Therapeutic implications with evidence levels and medications, diagnostic implications with clinical significance, molecular consequences and variant impact, per-variant oncogenicity prediction for simple variants, and ClinVar integration with star ratings
* **Ephemeral Oncogenicity Results**: Computed oncogenicity predictions remain visible only in the current page view and are cleared on a full page reload or in a new tab
* **Co-occurring Variants**: Flags CIViC molecular profiles (AND combinations of variants) when the patient has every component variant, shown as a collapsible box at the top of the Tx Implications cell
* **Pathway Diagrams**: Select variants and view them on the 10 TCGA oncogenic signaling pathways (Sanchez-Vega et al., *Cell* 2018), with oncogenes in pink, tumor suppressor genes in blue, the patient's variants outlined, and co-occurring variants connected
* **AI Synthesis**: a tumor-board brief per patient (reason for visit, case summary, genomic findings, cancer drivers, information gaps, ranked treatment options with trials, emerging evidence, references), with a linked table of contents. Pre-generated for the demo; see `public/data/synthesis/README.md` for how they're written
* **User Experience**: Responsive design optimized for clinical workflows, tooltips and help text for complex terminology, and an in-app How to Use guide

---

## Data Sources

* **Clinical data:** a FHIR R4 server (currently `https://hapi.fhir.org/baseR4`). Patients are found by MRN (`Patient.identifier` with system `http://elimu.io/fhir/mrn`) and loaded with `Patient/$everything`. The demo patients' mCODE transaction bundles, and the guide to creating new cases, are in `clinical-data/` (see `clinical-data/README.md`).
* **Genomic data:** the Genomics Operations reference server (`https://fhir-gen-ops.herokuapp.com`): subject operations for variants, structural variants, molecular consequences and Tx/Dx implications, with the patient's MRN as subject. For layering reported variants, the app also calls the server's `utilities` for the few reported variants that need it: `normalize-variant` (a reported variant's SPDI in each genome build) and `get-feature-coordinates` (a gene's or transcript's range), cached for the session.
* **Knowledge files** in `public/data/`: preset gene lists, the Dx crosswalk, the cancer-type codes for each tumor type, the CIViC co-occurrence knowledge base, the pathway diagrams, the OncoKB cancer gene list, and the pre-generated syntheses.

---

## Tech Stack

* **Framework**: Next.js 15 with React 19
* **Language**: TypeScript
* **Styling**: Tailwind CSS
* **UI Components**: Radix UI (tooltips), Lucide React (icons)
* **Data Processing**: PapaParse for CSV handling, Lodash utilities
* **Architecture**: Component-based with service layer abstraction

---

## Project Structure

```
src/
├── app/                    # Next.js app router
│   ├── globals.css         # Global styles
│   ├── howToUseContent.ts  # In-app How to Use guide
│   ├── layout.tsx          # Root layout
│   └── page.tsx            # Main page: banner, tabs, Genomics Explorer
├── components/             # Reusable UI components
│   ├── PatientBanner.tsx   # MRN, clinical status, Dx
│   ├── caseSummary/        # Case Summary sections, incl. ReportedVariantsSection.tsx
│   ├── SynthesisView.tsx   # AI Synthesis tab (Markdown rendering, table of contents)
│   ├── sidebar/            # Filter sidebar components
│   ├── ResultsTable.tsx    # Main results display
│   ├── resultsTableColumns.ts  # Table column information
│   ├── SearchForm.tsx      # Search input and controls
│   ├── SearchStatus.tsx    # Search progress and per-term status
│   ├── FeedbackForm.tsx    # User feedback collection
│   ├── EmailSubscription.tsx # Email signup
│   ├── pathways/           # Pathways view and SVG pathway diagram
│   ├── CooccurrenceBox.tsx # Co-occurring variants box in the Tx cell
│   ├── SelectionTray.tsx   # Selected-variants bar
│   └── *Cell.tsx           # Table cell components
├── services/               # Data layer services
│   ├── clinicalService.ts       # Clinical FHIR server (patients by MRN, $everything)
│   ├── cachedVariantService.ts  # Cached variant lookups
│   ├── variantService.ts        # Core variant API
│   ├── txService.ts             # Therapeutic implications
│   ├── dxService.ts             # Diagnostic implications
│   ├── mcService.ts             # Molecular consequences
│   ├── reportedVariantLookupService.ts  # normalize-variant / get-feature-coordinates lookups for reported variants
│   ├── cooccurrenceService.ts   # CIViC co-occurrence knowledge base and matching
│   ├── pathwayService.ts        # Pathway diagram loading
│   └── synthesisService.ts      # Loads a patient's synthesis
├── lib/                    # Shared logic: case summary, reported variants and layering, treatment timeline, display labels
├── types/                  # TypeScript type definitions
├── cancerFilter/           # Cancer type filtering logic
└── utils/                  # Utility functions

public/data/
├── MTB_KB_GeneLists.csv        # Cancer types and preset gene lists
├── MTB_TumorTypeCrosswalk.csv  # Recorded diagnosis (SNOMED CT) to tumor type (built by terminology/build_crosswalk.py)
├── MTB_TumorTypeImplicationCodes.csv  # Cancer-type codes (Disease Ontology, NCIt) that count as each tumor type
├── CIViC_Cooccurrence_KB.csv   # Co-occurrence profiles (built by scripts/build_civic_cooccurrence_kb.py)
├── OncoKB_CancerGeneList.tsv   # OncoKB cancer gene list (used when writing syntheses)
├── pathways/                   # 10 TCGA pathway diagrams (JSON)
└── synthesis/                  # Pre-generated syntheses <MRN>.md, search logs, and the README on how they're written

clinical-data/                  # Demo patients' mCODE bundles and build scripts; README on creating new cases

scripts/
└── build_civic_cooccurrence_kb.py  # Rebuilds CIViC_Cooccurrence_KB.csv from the CIViC API
```

---

## Usage

### Basic Search Workflow

1. **Choose a patient**: pick an MRN in the banner. The Dx is set from the patient's recorded primary cancer; override it with a generic tumor type if needed.
2. **Optionally Load a Preset Search List**: use the buttons under the search box to populate it:

   * `Actionable Gene List` and `Extended Gene List` for the Dx's tumor type
   * `Reported genes`: the genes of the patient's reported variants
3. **Enter Search Terms**:

   * Gene symbols: `BRAF, EGFR, TP53`
   * Genomic ranges: `NC_000007.14:55019016-55211628`
   * Mixed queries: `BRAF, NC_000007.14:55174721-55174820`
4. **Apply Filters**: Use `Filters` to open the sidebar and refine results with actionability (including `Reported variants`), molecular consequence, therapeutic implication, and diagnostic filters
5. **Review Results**: reported variants carry a Reported chip; their implications a Report tag

### AI Synthesis

Open the AI Synthesis tab and click `Generate synthesis`. For the demo, syntheses are pre-generated Markdown files (`public/data/synthesis/<MRN>.md`) written by the process in `public/data/synthesis/README.md`; `synthesisService.ts` could be replaced by a live AI call without changing the tab. Section headings get ids from their text, so the synthesis's table of contents links scroll within the tab.

### Co-occurring Variants

The Tx Implications cell starts with a **Potentially relevant co-occurring variants** box when the patient has every variant in a CIViC molecular profile that combines variants with AND (for example, `EGFR L858R AND EGFR T790M`). A component counts as present when the patient's therapeutic implications include either its CIViC variant ID or its single-variant CIViC molecular profile ID (Cat-VRS results are identified by molecular profile, e.g., `MET Amplification` = MP 266), so profiles involving copy number or other categorical variants are found only when Cat-VRS queries are enabled. Evidence is worded like other therapeutic implications; drugs are joined by CIViC's interaction type: combination `A + B`, substitutes `A or B`, sequential `A, then B`.

The profiles come from `public/data/CIViC_Cooccurrence_KB.csv`, built from the CIViC GraphQL API by `scripts/build_civic_cooccurrence_kb.py` (Python standard library only). Rerun it by hand to refresh the data, then commit the CSV:

```
python scripts/build_civic_cooccurrence_kb.py
```

It keeps AND-only profiles with predictive evidence (accepted or submitted) and skips profiles with OR, NOT/wildtype parts, or fusions.

### Pathways

1. Select variants with the checkboxes in the first column of the results table, or use `Select all shown`.
2. Click `View on pathways` in the bar at the bottom, or the `Pathways` tab.
3. Pathways are ranked by how many selected genes each diagram contains; the view opens on the top-ranked one.
4. Click a gene to see the patient's variants in it and jump back to the row with `Show in results table`.

The 10 pathway diagrams are stored as JSON in `public/data/pathways/` (see the README there), redrawn from Figure 2 of Sanchez-Vega et al., *Cell* 2018, with gene membership and oncogene / tumor suppressor roles from Table S3. Selections are kept only in the current page view.

### Oncogenicity Prediction Column

The `Oncogenicity Prediction` column is available for simple variants only.

1. Click `Compute prediction` on an SNV, MNV, or InDel row.
2. The app converts the displayed SPDI-style variant into HGVS and requests an oncogenicity prediction.
3. The table cell updates to a color-coded gauge based on the returned numeric score.
4. Click the gauge to open a detailed modal showing the overall score, overall interpretation, HGVS variant, original SPDI, evidence-line scoring, and caveats.
5. Use `View extended evidence details` in the modal to fetch the detailed evidence payload on demand.

Structural variants do not currently support oncogenicity prediction. Computed predictions are stored only in the current page view and are cleared on a full page reload or when the app is opened in a new tab.

### Search Examples

| Query Type     | Example                            | Description                 |
| -------------- | ---------------------------------- | --------------------------- |
| Gene Symbol    | `BRAF`                             | Find all BRAF variants      |
| Multiple Genes | `BRAF, EGFR, KRAS`                 | Search multiple genes       |
| Genomic Range  | `NC_000007.14:55019016-55211628`   | Search specific coordinates |
| Mixed Query    | `BRAF, NC_000017.11:7687550` | Combine different formats   |

### Filter System

**Preset Search Lists:**

* Populate the search box with cancer-type-specific predefined terms
* Terms may be gene symbols, genomic ranges, or a mix of both

**Sidebar Filters:**

* Actionability options: `Actionable, this tumor type`, `Actionable, any tumor type`, `Possibly actionable`, `Reported variants`, or `None`
* Evidence levels, medications, implications
* Molecular consequences and impact levels
* Diagnostic significance and ClinVar ratings
* Actionable, this tumor type matches the implication's cancer-type code against the tumor type's codes

---

## Development

### Key Components

* **PatientBanner**: MRN selection, clinical status and Dx
* **CaseSummary / ReportedVariantsSection**: The clinical record and the lab's reported variants
* **SynthesisView**: Renders a patient's synthesis
* **SearchForm**: Handles user input, preset search-term loading, and search execution
* **ResultsTable**: Displays variants with expandable details
* **FilterSidebar**: Advanced filtering interface, including actionability controls
* **CooccurrenceBox**: Collapsible co-occurring variants box at the top of the Tx Implications cell
* **SelectionTray**: Bar showing the selected variants, with `View on pathways`
* **PathwaysView / PathwayDiagram**: Ranked pathway list, selected variants, SVG diagram, and gene details

### Services Architecture

The application uses a service layer pattern:

* **clinicalService**: Clinical FHIR server access
* **variantService**: Core API communication
* **reportedVariantLookupService**: Genomics server utility lookups for reported variants
* **synthesisService**: Loads pre-generated syntheses
* **cachedVariantService**: Caching layer for performance
* **txService/dxService/mcService**: Specialized implication handlers
* **cacheService**: Generic caching utilities
* **cooccurrenceService**: Loads the CIViC co-occurrence knowledge base and matches it against the patient's results
* **pathwayService**: Loads the stored pathway diagrams

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
