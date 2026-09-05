# Molecular Tumor Board (MTB) — Variant & Implication Explorer

Live demo: [https://moleculartb.netlify.app/](https://moleculartb.netlify.app/)

A web application that helps clinicians and researchers quickly review **therapeutic (Tx)**, **diagnostic (Dx)**, and **molecular consequence** data for cancer variants. The platform supports filtering by cancer type, actionability levels, and phenotype to streamline Molecular Tumor Board (MTB) review workflows.

![MTB App Interface](public/elimu.png)

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
* **Functionality**: The application allows read-only searches of genetic variants and their annotations. Clinicians can enter ranges (gene names or genomic coordinates), select cancer types, and apply actionability filters to quickly narrow down relevant results.
* **Outcome**: Results are displayed in a structured table, showing molecular consequences, diagnostic significance, and therapeutic implications, with links to external reference databases (ClinVar, CIViC, SnpEff).

This workflow ensures that clinicians can:

* Review variants relevant to specific cancer types
* Filter by actionability (evidence level A, tumor-specific or other tumor types)
* Analyze therapeutic implications and drug responses
* Assess diagnostic significance and molecular consequences
* Rely on evidence-linked annotations during MTB discussions

---

## Core Features

* **Intelligent Search**: Search by gene symbols (e.g., `BRAF`, `EGFR`) or genomic ranges (e.g., `NC_000007.14:55019016-55211628`)
* **Advanced Filtering**: Cancer type selection, actionability levels, evidence levels, medications, implications, and molecular consequences
* **Real-time Results**: Incremental loading with search status indicators and cancellation support
* **Comprehensive Data Display**: Therapeutic implications with evidence levels and medications, diagnostic implications with clinical significance, molecular consequences and variant impact, ClinVar integration with star ratings
* **User Experience**: Responsive design optimized for clinical workflows, expandable result rows with detailed information, tooltips and help text for complex terminology

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
├── app/                     # Next.js app router
│   ├── globals.css         # Global styles
│   ├── layout.tsx          # Root layout
│   └── page.tsx            # Main application page
├── components/             # Reusable UI components
│   ├── sidebar/            # Filter sidebar components
│   ├── ResultsTable.tsx    # Main results display
│   ├── SearchForm.tsx      # Search input and controls
│   ├── FeedbackForm.tsx    # User feedback collection
│   ├── EmailSubscription.tsx # Email signup
│   └── *Cell.tsx           # Table cell components
├── services/               # Data layer services
│   ├── cachedVariantService.ts  # Cached variant lookups
│   ├── variantService.ts        # Core variant API
│   ├── txService.ts            # Therapeutic implications
│   ├── dxService.ts            # Diagnostic implications
│   └── mcService.ts            # Molecular consequences
├── types/                  # TypeScript type definitions
│   └── variants.ts         # Variant and implication types
├── cancerFilter/           # Cancer type filtering logic
├── utils/                  # Utility functions
└── public/                # Static assets
```

---

## Usage

### Basic Search Workflow

1. **Select Cancer Type**: Choose from the dropdown to filter relevant variants
2. **Set Actionability**: Use radio buttons to focus on:

   * "Actionable, this tumor type" (A-level evidence for selected cancer)
   * "Actionable, other tumor type" (A-level evidence for other cancers)
   * "Possibly actionable" (non-A level evidence)
3. **Enter Search Terms**:

   * Gene symbols: `BRAF, EGFR, TP53`
   * Genomic ranges: `NC_000007.14:55019016-55211628`
   * Mixed queries: `BRAF V600E, NC_000007.14:55174721-55174820`
4. **Apply Advanced Filters**: Use the filter sidebar for granular control
5. **Review Results**: Expand rows to see detailed implications

### Search Examples

| Query Type     | Example                            | Description                 |
| -------------- | ---------------------------------- | --------------------------- |
| Gene Symbol    | `BRAF`                             | Find all BRAF variants      |
| Multiple Genes | `BRAF, EGFR, KRAS`                 | Search multiple genes       |
| Genomic Range  | `NC_000007.14:55019016-55211628`   | Search specific coordinates |
| Mixed Query    | `BRAF V600E, NC_000017.11:7687550` | Combine different formats   |

### Filter System

**Quick Filters (Radio Buttons):**

* Pre-configured evidence level and phenotype matching
* Automatically applied based on cancer type selection

**Advanced Filters (Sidebar):**

* Evidence levels, medications, implications
* Molecular consequences and impact levels
* Diagnostic significance and ClinVar ratings
* Custom phenotype matching

---

## Development

### Key Components

* **SearchForm**: Handles user input and search execution
* **ResultsTable**: Displays variants with expandable details
* **FilterSidebar**: Advanced filtering interface
* **CancerSelect**: Cancer type selection with phenotype loading
* **ActionableCheckBoxes**: Quick actionability filtering

### Services Architecture

The application uses a service layer pattern:

* **variantService**: Core API communication
* **cachedVariantService**: Caching layer for performance
* **txService/dxService/mcService**: Specialized implication handlers
* **cacheService**: Generic caching utilities

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
