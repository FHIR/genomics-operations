export interface HowToUseLinkItem {
    textBefore?: string;
    linkLabel: string;
    linkUrl: string;
    textAfter?: string;
}

export interface HowToUseSection {
    title: string;
    items: Array<string | HowToUseLinkItem>;
}

// Top-level groups match the tabs below the banner; a group without a title holds what applies to all of them
export interface HowToUseGroup {
    title?: string;
    sections: HowToUseSection[];
}

export const HOW_TO_USE_TITLE = 'How to use the Molecular Tumor Board app';

export const HOW_TO_USE_GROUPS: HowToUseGroup[] = [
    {
        sections: [
        {
            title: 'Introduction',
            items: [
                'This is a proof-of-concept app, designed as a general purpose genetic data viewer, customized for use by a Molecular Tumor Board.',
                'The primary objective of this open-source app is enable rapid prototyping and experimentation with advanced FHIR capabilities needed to support MTB functionality.',
                {
                    textBefore: 'Source code for the app is: ',
                    linkLabel: 'here',
                    linkUrl: 'https://github.com/FHIR/genomics-operations/tree/main/genomics-apps/mtb-app',
                    textAfter: '.'
                },
            ],
        },
        {
            title: 'Patient Banner',
            items: [
                'The banner at the top of the page stays in view as you scroll. Choose a patient from its MRN dropdown, or pick Enter custom MRN at the bottom of the list.',
                'Below the banner, the Case Summary, Genomic Findings and Synthesis tabs switch views. The app starts with no patient selected: choose one from the MRN box in the banner. The first patient you choose opens on Case Summary if they have clinical data, otherwise on Genomic Findings; after that your tab stays as you switch patients. Searches and results are kept while you look at the other tabs.',
                'Patients are grouped into With clinical data (patients the app found on the clinical FHIR server) and Genomic data only. Switching patients clears the previous patient\'s genomic results; the search terms stay so you can rerun them.',
                'For a patient with clinical data, the banner shows name, age, sex and date of birth, plus the latest stage (with the year it was assigned), ECOG performance status and cancer disease status from the patient\'s record. Hover over an item for its date and details.',
                'When the patient is enrolled in a research study, an outlined On study chip with a flask icon follows the disease status, showing the study\'s ClinicalTrials.gov number (NCT…). Hover over it for the study title and enrollment date; click it to open the study on ClinicalTrials.gov. A patient in more than one study shows On 2 studies, with each study listed on hover. Screening and past enrollments aren\'t shown.',
                'Stage is shown as recorded and is never recalculated. When metastases were recorded after a staging, the Case Summary notes it on that stage, so an early stage isn\'t mistaken for the current extent of disease.',
                'The Dx dropdown starts with the patient\'s own primary cancer, as recorded. You can override it with a generic tumor type (Breast Cancer, Hepatocellular Carcinoma, Non-Small Cell Lung Cancer). A patient without clinical data starts with no Dx.',
                'The Dx sets the cancer-specific gene-list presets, the Actionable, this tumor type filter, and the tumor type sent to the oncogenicity predictor. The predictor gets no tumor type when it has no matching one, for example for Hepatocellular Carcinoma.',
                {
                    textBefore: 'A recorded diagnosis maps to a tumor type when its SNOMED CT code is in that tumor type\'s ',
                    linkLabel: 'VSAC',
                    linkUrl: 'https://vsac.nlm.nih.gov/',
                    textAfter: ' value set (exact code match). Other diagnoses still appear in the Dx dropdown but set no presets or predictor tumor type.',
                },
            ],
        },
        ],
    },
    {
        title: 'Case Summary',
        sections: [
        {
            title: 'Overview',
            items: [
                'The Case Summary shows the patient\'s clinical record from the FHIR server: cancer diagnosis and staging, a treatment timeline, lab and pathology results, imaging observations, assessments, vitals and social history, other conditions, allergies, family history and reports. Administrative records (practitioners, encounters, specimens and the like) aren\'t shown.',
                'The bar of section links at the top of the Case Summary stays in view as you scroll; click one to jump to that section. The section you are viewing is highlighted.',
                'Labels use the text the record\'s author wrote, else the code\'s display name, else the code itself.',
                'Longer lists in the smaller cards scroll within their card, so the page layout stays the same however much data a patient has.',
            ],
        },
        {
            title: 'Cancer Diagnosis & Staging',
            items: [
                'Shows the primary cancer as recorded, with its newest stage, histology, body site, onset and status, then each stage group with its T, N and M categories, then any metastases.',
                'Stage is shown as recorded and is never recalculated. When metastases were recorded after a staging, that stage carries a note listing them.',
            ],
        },
        {
            title: 'Treatment Timeline',
            items: [
                'Shows every medication (active first, then inactive, each alphabetical), radiotherapy courses, procedures, cancer disease status and the dates of diagnosis and metastases, from the earliest event to today. Hover over a bar or marker for its dates and details.',
                'Medication bars span the recorded dosing period; active ones run to the dashed Today line. When an order has no dosing period, its bar starts at the order date and its label says so. Disease status uses the record\'s own terms (Improved, Stable, Worsened), not RECIST response categories.',
            ],
        },
        {
            title: 'Lab & Pathology Results',
            items: [
                'Everything recorded as a laboratory result, including tumor markers and pathology findings, one row per test with the latest value, reference range and flag. Sort by test name or date, and click a row for its full history; tests with 3 or more numeric values get a trend line.',
            ],
        },
        ],
    },
    {
        title: 'Genomic Findings',
        sections: [
        {
            title: 'Search Basics',
            items: [
                'Choose a patient in the banner before searching.',
                'Enter one or more gene symbols separated by commas, for example: BRAF, EGFR, ALK.',
                'You can also search by genomic range using zero-based RefSeq:start-end format, for example: NC_000007.14:55019016-55211628.',
                'When the Dx has gene lists, use the Actionable Gene List or Extended Gene List preset buttons under the search box to load them. Hepatocellular Carcinoma has only an Extended Gene List.',
            ],
        },
        {
            title: 'Search Modes',
            items: [
                'The app retrieves both simple variants and structural variants for the selected patient and searched region, along with select annotations, including therapeutic options.',
                {
                    textBefore: 'When ',
                    linkLabel: 'Cat-VRS',
                    linkUrl: 'https://cat-vrs.ga4gh.org/en/stable/',
                    textAfter: ' queries is checked, therapeutic implication queries include constraint checking against Cat-VRS encoded knowledge.',
                },
            ],
        },
        {
            title: 'Understanding the Results Table',
            items: [
                'Results are grouped by the searched range or gene, each group under a band naming it and its variant count, and are sorted by Range by default. Click the arrow beside a column name to sort by that column (up for A-Z, down for Z-A).',
                'Simple variants are labeled as SNV, MNV, or InDel.',
                'Structural variants display the DNA change type first, followed by the genomic location in parentheses, and copy number when available.',
                'The Variant column may also show the first protein change below the main variant label when a molecular consequence is available.',
            ],
        },
        {
            title: 'Filtering and Table Controls',
            items: [
                'Use Filters to open the sidebar filters for actionability, molecular consequences, therapeutic implications, and diagnostic implications.',
                'Use Customize table to show or hide columns, change their order, resize them, and apply per-column filtering and sorting.',
            ],
        },
        {
            title: 'Co-occurring Variants',
            items: [
                'When the patient has every variant named in a CIViC molecular profile that combines variants (for example, EGFR L858R AND EGFR T790M), the Tx Implications cell starts with a Potentially relevant co-occurring variants box. Click it to open or close it.',
                'Each profile name links to CIViC, followed by its evidence in the same format as other therapeutic implications. Drugs given together are joined with +, interchangeable drugs with or, and drugs given in sequence with then.',
                'Matching uses the CIViC variant and molecular profile identifiers returned with each variant\'s therapeutic implications. Profiles involving copy number changes or other categorical variants are only found when Cat-VRS queries is checked. Profiles that use OR, NOT or wildtype, and fusion profiles, are not checked yet.',
                'The box is not affected by the Filters options.',
            ],
        },
        {
            title: 'Selecting Variants and Pathways',
            items: [
                'Use the checkboxes in the first column of the results table, or Select all shown, to select variants. Selected variants stay selected when you change filters; the bar at the bottom shows how many are hidden by current filters.',
                'Click View on pathways (or the Pathways tab) to see the selected variants on the 10 TCGA oncogenic signaling pathways, ranked by how many of your selected genes each diagram contains.',
                'In the diagrams, pink boxes are oncogenes and blue boxes are tumor suppressor genes. Boxes with a thick dark outline contain your selected variants. A dashed purple arc connects genes in a matching co-occurrence profile.',
                'Click a gene to see the patient\'s variants in it, their therapeutic implications, and a link back to the row in the results table.',
                {
                    textBefore: 'The diagrams are redrawn from Figure 2 of ',
                    linkLabel: 'Sanchez-Vega et al., Cell 2018',
                    linkUrl: 'https://doi.org/10.1016/j.cell.2018.03.035',
                    textAfter: ', with gene roles from the paper\'s Table S3. Selections are kept only in the current page view.',
                },
            ],
        },
        {
            title: 'Oncogenicity Prediction Column',
            items: [
                'The Oncogenicity Prediction column applies only to simple variants (SNVs, MNVs, and InDels). Structural variants currently show that oncogenicity prediction is not available.',
                'For a simple variant, click Compute prediction to run the oncogenicity predictor. The cell then shows a color-coded gauge, and selecting the gauge opens a detailed modal.',
                'The first prediction after a period of inactivity can take up to a minute while the prediction server wakes up. If a prediction fails, it is retried automatically; if it still fails, the cell shows Prediction failed (hover for the reason) with a Retry button.',
                'The detailed modal shows the overall score, overall prediction, submitted variant, evidence-line scores, and caveats. You can also request extended evidence details from the same modal.',
                'Computed oncogenicity predictions are kept only in the current page view. They are cleared if you reload the page or open the app in a new tab, and when the Dx changes, since the prediction depends on the tumor type.',
            ],
        },
        {
            title: 'Links and Evidence',
            items: [
                'Therapeutic implications can include external links such as CIViC entries and, for clinical trial matches, a Clinical Trial link.',
                'Diagnostic and therapeutic evidence levels are shown in brackets within the table so you can compare the strength of support quickly.',
                'If no implication or consequence is available for a variant, the table shows <none found>.',
            ],
        },
        {
            title: 'Practical Notes',
            items: [
                'Large multi-gene searches are processed one search term at a time. The Status row under the search box shows a summary, a progress bar and a chip for each term; for long searches the chips stop after about two rows, and Show all lists every term in a scrolling box.',
                'You can cancel an in-progress search with Cancel search in the Status row. Hide collapses the search box to one line showing the search terms and a status summary.',
            ],
        },
        {
            title: 'Example',
            items: [
                'Patient Larry Lung (MRN L2345) is a 66 yo male, never-smoker, who presented May 10, 2023 with persistent intermittent dry cough, found to have a RLL lung mass on CXR. Chest CT and PET showed the RLL lung mass with right hilar adenopathy and additional pulmonary and bony metastases. Brain MRI with no evidence of intracranial metastases. Fine needle aspiration of the RLL lung mass showed adenocarcinoma, as did the left iliac biopsy specimen. Tumor genomic profiling has been performed and the clinical team now wants to see if the patient has any actionable variants.',
                'Select MRN "L2345" in the banner and open the Genomic Findings tab. The banner shows his stage, ECOG and latest disease status, and the Dx is set from his record to Adenocarcinoma of right lower lobe of lung, which maps to Non-Small Cell Lung Cancer.',
                'Click Actionable Gene List to populate the search box with the predefined NSCLC actionable gene list, then run the search.',
                'Use the Filters button to confirm that Actionable, this tumor type is selected by default, or switch actionability to Actionable, any tumor type, Possibly actionable, or None.',
                'Use the same sidebar to apply additional filters for molecular consequences, therapeutic implications, and diagnostic implications.',
                'With Cat-VRS queries checked, the Tx Implications for EGFR L858R start with Potentially relevant co-occurring variants, including EGFR L858R AND MET Amplification.',
                'Tick the EGFR and MET variants and click View on pathways to see them on the RTK/RAS pathway, joined by a co-occurring variants arc.',
            ],
        },
        ],
    },
    {
        title: 'Synthesis',
        sections: [
            {
                title: 'Generating a Synthesis',
                items: [
                    'Open the Synthesis tab and click Generate synthesis. After a short progress display, the synthesis appears: a case summary, genomic findings, information gaps and treatment options, including recruiting clinical trials, with references.',
                    'In this demo the syntheses are pre-generated for the four example patients (L2345, HCC1143, MRN1234 and TCGA-DD-A1EH) and labeled as such; other patients show that no synthesis is available. A future version could generate them live with an AI model.',
                    'A generated synthesis stays on the tab while you look at other tabs, and is cleared when you choose a different patient.',
                    'The synthesis is a research demonstration based on synthetic patients and is not for clinical use.',
                ],
            },
        ],
    },
];
