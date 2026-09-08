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
            'The primary objective of this open-source app is to enable experimentation with advanced FHIR capabilities.',
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
            'If you select a Cancer Type, you can then select one of the "Filter by actionability" options, to load a predefined gene list.',
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
                textAfter: ' queries (experimental) is checked, therapeutic implication queries include constraint checking against Cat-VRS encoded knowledge.',
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
            'Use Filter Results to open the implication filters for molecular consequences, therapeutic implications, and diagnostic implications.',
            'Use Customize Table to show or hide columns, change their order, resize them, and apply per-column filtering and sorting.',
            'Some columns are hidden by default, such as Oncogenicity Prediction.',
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
            'To limit search to just actionable genes, and limit results to only those with high grade evidence specific for this cancer type, select "Actionable, this tumor type". ',
            'All variants in specified ranges are returned, but filtered based on actionability button. Unclick the actionability button to remove all filters.',
            'Use the "Filter Results" button to apply additional filters for molecular consequences, therapeutic implications, and diagnostic implications.',
        ],
    },

];
