# Molecular Tumor Board (MTB) — Variant & Implication Explorer

Live demo: [https://mtb-app-elimu1.vercel.app/](https://mtb-app-elimu1.vercel.app/)

A web application that helps clinicians and researchers quickly review **therapeutic (Tx)**, **diagnostic (Dx)**, and **molecular consequence** data for cancer variants. The platform supports cancer-type-specific search presets, actionability filtering, and phenotype-aware review to streamline Molecular Tumor Board (MTB) workflows.

---

## Table of Contents

* [Overview](#overview)
* [Core Features](#core-features)
* [Tech Stack](#tech-stack)
* [Project Structure](#project-structure)
* [Usage](#usage)
* [Development](#development)
* [License](#license)

---

## Overview

The **MTB Variant & Implication Explorer** was designed to support oncologists during Molecular Tumor Board meetings by providing a fast and reliable way to explore genomic profiling results. The tool was built around a typical MTB workflow:

* **Scenario**: A patient’s tumor genomic profiling results are presented at the MTB. Clinicians need to review variants, assess their therapeutic and diagnostic implications, and decide on treatment recommendations.
* **Functionality**: The application allows read-only searches of genetic variants and their annotations. Clinicians can enter gene symbols or genomic ranges, load predefined gene lists for a selected cancer type, and apply actionability filters from the results sidebar.
* **Outcome**: Results are displayed in a structured table, showing molecular consequences, diagnostic significance, and therapeutic implications, with links to external reference databases (e.g. ClinVar, CIViC, SnpEff).

This workflow ensures that clinicians can:

* Review variants relevant to specific cancer types
* Load predefined gene lists for a selected cancer type
* Filter by actionability, evidence level, and phenotype context
* Analyze therapeutic implications and drug responses
* Assess diagnostic significance and molecular consequences
* Rely on evidence-linked annotations during MTB discussions

---

## Core Features

* **Intelligent Search**: Search by gene symbols (e.g., `BRAF`, `EGFR`) or genomic ranges (e.g., `NC_000007.14:55019016-55211628`)
* **Guided Search Presets**: Cancer-type-specific Actionable Gene List and Extended Gene List buttons populate the search box with predefined terms
* **Advanced Filtering**: Sidebar actionability, evidence levels, medications, implications, and molecular consequences
* **Real-time Results**: Incremental loading with search status indicators and cancellation support
* **Comprehensive Data Display**: Therapeutic implications with evidence levels and medications, diagnostic implications with clinical significance, molecular consequences and variant impact, per-variant oncogenicity prediction for simple variants, and ClinVar integration with star ratings
* **Ephemeral Oncogenicity Results**: Computed oncogenicity predictions remain visible only in the current page view and are cleared on a full page reload or in a new tab
* **User Experience**: Responsive design optimized for clinical workflows, expandable result rows with detailed information, tooltips and help text for complex terminology
* **Co-occurring Variants**: Flags CIViC molecular profiles (AND combinations of variants) when the patient has every component variant, shown as a collapsible box at the top of the Tx Implications cell
* **Pathway Diagrams**: Select variants and view them on the 10 TCGA oncogenic signaling pathways (Sanchez-Vega et al., *Cell* 2018), with oncogenes in pink, tumor suppressor genes in blue, the patient's variants outlined, and co-occurring variants connected

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
│   ├── howToUseContent.ts  # Overview of app use
│   ├── layout.tsx          # Root layout
│   └── page.tsx            # Main application page
├── components/             # Reusable UI components
│   ├── sidebar/            # Filter sidebar components
│   ├── ResultsTable.tsx    # Main results display
│   ├── resultsTableColumns.ts  # Table column information
│   ├── SearchForm.tsx      # Search input and controls
│   ├── FeedbackForm.tsx    # User feedback collection
│   ├── EmailSubscription.tsx # Email signup
│   ├── pathways/           # Pathways view and SVG pathway diagram
│   ├── CooccurrenceBox.tsx # Co-occurring variants box in the Tx cell
│   ├── SelectionTray.tsx   # Selected-variants bar
│   └── *Cell.tsx           # Table cell components
├── services/               # Data layer services
│   ├── cachedVariantService.ts  # Cached variant lookups
│   ├── variantService.ts        # Core variant API
│   ├── txService.ts            # Therapeutic implications
│   ├── dxService.ts            # Diagnostic implications
│   ├── mcService.ts            # Molecular consequences
│   ├── cooccurrenceService.ts  # CIViC co-occurrence knowledge base and matching
│   └── pathwayService.ts       # Pathway diagram loading
├── lib/                    # Shared helpers (e.g., variant display labels)
├── types/                  # TypeScript type definitions
│   └── variants.ts         # Variant and implication types
├── cancerFilter/           # Cancer type filtering logic
├── utils/                  # Utility functions
└── public/                # Static assets

public/data/
├── MTB_KB_GeneLists.csv        # Cancer types and preset gene lists
├── CIViC_Cooccurrence_KB.csv   # Co-occurrence profiles (built by scripts/build_civic_cooccurrence_kb.py)
└── pathways/                   # 10 TCGA pathway diagrams (JSON)

scripts/
└── build_civic_cooccurrence_kb.py  # Rebuilds CIViC_Cooccurrence_KB.csv from the CIViC API
```

---

## Usage

### Basic Search Workflow

1. **Select Cancer Type**: Choose from the dropdown to enable cancer-specific search presets and phenotype-aware filtering
2. **Optionally Load a Preset Search List**: Use the gray search panel buttons to populate the search box for the selected cancer type:

   * `Actionable Genes`
   * `Extended Gene List`
3. **Enter Search Terms**:

   * Gene symbols: `BRAF, EGFR, TP53`
   * Genomic ranges: `NC_000007.14:55019016-55211628`
   * Mixed queries: `BRAF, NC_000007.14:55174721-55174820`
4. **Apply Filters**: Use `Filter Results` to open the sidebar and refine results with actionability, molecular consequence, therapeutic implication, and diagnostic filters
5. **Review Results**: Expand rows to see detailed implications

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

* Actionability options: `Actionable, this tumor type`, `Actionable, any tumor type`, `Possibly actionable`, or `None`
* Evidence levels, medications, implications
* Molecular consequences and impact levels
* Diagnostic significance and ClinVar ratings
* Custom phenotype matching

---

## Development

### Key Components

* **SearchForm**: Handles user input, preset search-term loading, and search execution
* **ResultsTable**: Displays variants with expandable details
* **FilterSidebar**: Advanced filtering interface, including actionability controls
* **CancerSelect**: Cancer type selection with phenotype loading
* **CooccurrenceBox**: Collapsible co-occurring variants box at the top of the Tx Implications cell
* **SelectionTray**: Bar showing the selected variants, with `View on pathways`
* **PathwaysView / PathwayDiagram**: Ranked pathway list, selected variants, SVG diagram, and gene details

### Services Architecture

The application uses a service layer pattern:

* **variantService**: Core API communication
* **cachedVariantService**: Caching layer for performance
* **txService/dxService/mcService**: Specialized implication handlers
* **cacheService**: Generic caching utilities
* **cooccurrenceService**: Loads the CIViC co-occurrence knowledge base and matches it against the patient's results
* **pathwayService**: Loads the stored pathway diagrams

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
