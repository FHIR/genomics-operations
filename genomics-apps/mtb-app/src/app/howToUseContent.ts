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

export const HOW_TO_USE_TITLE = 'How to use the Molecular Tumor Board Genetic Data Viewer';

export const HOW_TO_USE_SECTIONS: HowToUseSection[] = [
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
        title: 'Search Basics',
        items: [
            'Choose an MRN before searching. The default MRN is useful for demos, but you can select a different patient or enter a custom MRN.',
            'Enter one or more gene symbols separated by commas, for example: BRAF, EGFR, ALK.',
            'You can also search by genomic range using zero-based RefSeq:start-end format, for example: NC_000007.14:55019016-55211628.',
            'After you select a Cancer Type, you can use the Actionable Genes or Extended Gene List buttons in the search panel to load predefined gene lists.',
        ],
    },
    {
        title: 'Search Modes',
        items: [
            'The app retrieves both simple variants and structural variants for the selected patient and searched region, along with select annotations, including therapeutic options.',
            {
                textBefore: 'When Enable ',
                linkLabel: 'Cat-VRS',
                linkUrl: 'https://cat-vrs.ga4gh.org/en/stable/',
                textAfter: ' queries is checked, therapeutic implication queries include constraint checking against Cat-VRS encoded knowledge.',
            },
        ],
    },
    {
        title: 'Understanding the Results Table',
        items: [
            'Results are grouped by the searched range or gene and are sorted by Range by default.',
            'Simple variants are labeled as SNV, MNV, or InDel.',
            'Structural variants display the DNA change type first, followed by the genomic location in parentheses, and copy number when available.',
            'The Variant column may also show the first protein change below the main variant label when a molecular consequence is available.',
        ],
    },
    {
        title: 'Filtering and Table Controls',
        items: [
            'Use Filter Results to open the sidebar filters for actionability, molecular consequences, therapeutic implications, and diagnostic implications.',
            'Use Customize Table to show or hide columns, change their order, resize them, and apply per-column filtering and sorting.',
        ],
    },
    {
        title: 'Oncogenicity Prediction Column',
        items: [
            'The Oncogenicity Prediction column applies only to simple variants (SNVs, MNVs, and InDels). Structural variants currently show that oncogenicity prediction is not available.',
            'For a simple variant, click Compute prediction to run the oncogenicity predictor. The cell then shows a color-coded gauge, and selecting the gauge opens a detailed modal.',
            'The detailed modal shows the overall score, overall prediction, HGVS conversion, original SPDI, evidence-line scores, and caveats. You can also request extended evidence details from the same modal.',
            'Computed oncogenicity predictions are kept only in the current page view. They are cleared if you reload the page or open the app in a new tab.',
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
            'Large multi-gene searches are processed one search term at a time, and Search Status shows progress for each term.',
            'You can cancel an in-progress search from the Search Status panel.',
        ],
    },
    {
        title: 'Example',
        items: [
            'Patient Larry Lung (MRN L2345) is a 66 yo male, never-smoker, who presented May 10, 2023 with persistent intermittent dry cough, found to have a RLL lung mass on CXR. Chest CT and PET showed the RLL lung mass with right hilar adenopathy and additional pulmonary and bony metastases. Brain MRI with no evidence of intracranial metastases. Fine needle aspiration of the RLL lung mass showed adenocarcinoma, as did the left iliac biopsy specimen. Tumor genomic profiling has been performed and the clinical team now wants to see if the patient has any actionable variants.',
            'Select MRN "L2345" and Cancer Type "NSCLC".',
            'Click Actionable Genes to populate the search box with the predefined NSCLC actionable gene list, then run the search.',
            'Use the Filter Results button to confirm that Actionable, this tumor type is selected by default, or switch actionability to Actionable, any tumor type, Possibly actionable, or None.',
            'Use the same sidebar to apply additional filters for molecular consequences, therapeutic implications, and diagnostic implications.',
        ],
    },

];
