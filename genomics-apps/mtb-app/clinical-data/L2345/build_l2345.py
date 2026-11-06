"""Generate the mCODE transaction bundle for Larry Lung (MRN L2345).

Usage: python build_l2345.py [output_path]
Writes L2345-transaction-bundle.json next to this script by default. The clinical
story behind the data is in L2345-scenario.md; edit this script, not the JSON.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # shared mcode_helpers.py

from mcode_helpers import (  # noqa: E402
    ICD10, LEFT, LOINC, MCODE, RIGHT, RXN, SCT, USCORE, add, cc, condition_base, ext,
    laterality, obs_cat, qty, specimen_type, write_bundle,
)

OUT = sys.argv[1] if len(sys.argv) > 1 else Path(__file__).with_name('L2345-transaction-bundle.json')

PID = 'mtb-l2345'
PAT = {'reference': f'Patient/{PID}', 'display': 'Larry Lung'}
ONC = {'reference': 'Practitioner/mtb-l2345-oncologist', 'display': 'Dr. Olivia Oncologist'}
PATH = {'reference': 'Practitioner/mtb-l2345-pathologist', 'display': 'Dr. Paul Pathologist'}
RAD = {'reference': 'Practitioner/mtb-l2345-radiologist', 'display': 'Dr. Rita Radiologist'}
RADONC = {'reference': 'Practitioner/mtb-l2345-radiation-oncologist', 'display': 'Dr. Ray Radonc'}
PRIMARY = {'reference': 'Condition/mtb-l2345-primary-cancer', 'display': 'Adenocarcinoma of right lower lobe of lung'}


# ---------------------------------------------------------------- Patient & practitioners

add({
    'resourceType': 'Patient', 'id': PID,
    'extension': [{'url': 'http://hl7.org/fhir/us/core/StructureDefinition/us-core-birthsex', 'valueCode': 'M'}],
    'identifier': [{
        'use': 'usual',
        'type': cc('http://terminology.hl7.org/CodeSystem/v2-0203', 'MR', 'Medical Record Number'),
        'system': 'http://elimu.io/fhir/mrn',
        'value': 'L2345',
    }],
    'active': True,
    'name': [{'use': 'official', 'family': 'Lung', 'given': ['Larry']}],
    'gender': 'male',
    'birthDate': '1957-03-14',
}, [MCODE + 'mcode-cancer-patient'])

for pid, family, given, specialty_code, specialty in [
    ('oncologist', 'Oncologist', 'Olivia', '394593009', 'Medical oncology'),
    ('pathologist', 'Pathologist', 'Paul', '394915009', 'General pathology (specialty)'),
    ('radiologist', 'Radiologist', 'Rita', '394914008', 'Radiology - specialty'),
    ('radiation-oncologist', 'Radonc', 'Ray', '419815003', 'Radiation oncology'),
]:
    add({
        'resourceType': 'Practitioner', 'id': f'mtb-l2345-{pid}',
        'identifier': [{'system': 'http://elimu.io/fhir/practitioner-id', 'value': f'mtb-{pid}'}],
        'name': [{'family': family, 'given': [given], 'prefix': ['Dr.']}],
        'qualification': [{'code': cc(SCT, specialty_code, specialty)}],
    }, [USCORE + 'us-core-practitioner'])

# ---------------------------------------------------------------- Specimens

SPEC_LUNG = add({
    'resourceType': 'Specimen', 'id': 'mtb-l2345-specimen-rll-fna',
    'status': 'unavailable',
    'type': specimen_type('TISS', 'Tissue', '122610009', 'Specimen from lung obtained by biopsy'),
    'note': [{'text': 'RLL mass, CT-guided FNA (cytology smears and cell block)'}],
    'subject': PAT,
    'collection': {
        'collector': RAD, 'collectedDateTime': '2023-05-25',
        'method': cc(SCT, '432231006', 'Fine needle aspiration biopsy of lung'),
        'bodySite': cc(SCT, '266005', 'Structure of lower lobe of right lung'),
    },
}, [MCODE + 'mcode-human-specimen'])

SPEC_ILIUM = add({
    'resourceType': 'Specimen', 'id': 'mtb-l2345-specimen-left-iliac-core',
    'status': 'unavailable',
    'type': specimen_type('BON', 'Bone', '16213011000119107', 'Specimen from bone obtained by core needle biopsy'),
    'note': [{'text': 'Left iliac bone lesion, CT-guided core biopsy'}],
    'subject': PAT,
    'collection': {
        'collector': RAD, 'collectedDateTime': '2023-05-30',
        'method': cc(SCT, '432371004', 'CT guided biopsy of bone'),
        'bodySite': cc(SCT, '773840006', 'Bone structure of left ilium'),
    },
}, [MCODE + 'mcode-human-specimen'])

SPEC_LIVER = add({
    'resourceType': 'Specimen', 'id': 'mtb-l2345-specimen-liver-core',
    'status': 'available',
    'type': specimen_type('TISS', 'Tissue', '258415003', 'Biopsy specimen'),
    'note': [{'text': 'Liver segment 6 lesion, CT-guided core biopsy'}],
    'subject': PAT,
    'collection': {
        'collector': RAD, 'collectedDateTime': '2026-10-06',
        'method': cc(SCT, '407971000119109', 'CT guided percutaneous needle biopsy of liver'),
        'bodySite': cc(SCT, '10200004', 'Liver structure'),
    },
}, [MCODE + 'mcode-human-specimen'])

# ---------------------------------------------------------------- Diagnosis

STAGE_GROUP = {'reference': 'Observation/mtb-l2345-tnm-clinical-stage-group'}

add({
    'resourceType': 'Condition', 'id': 'mtb-l2345-primary-cancer',
    'extension': [ext('mcode-histology-morphology-behavior', 'valueCodeableConcept',
                      cc(SCT, '1187332001', 'Adenocarcinoma'))],
    **condition_base(),
    'code': {'coding': [
        {'system': SCT, 'code': '254626006', 'display': 'Adenocarcinoma of lung'},
        {'system': ICD10, 'code': 'C34.31', 'display': 'Malignant neoplasm of lower lobe, right bronchus or lung'},
    ], 'text': 'Adenocarcinoma of right lower lobe of lung'},
    'bodySite': [{
        'extension': laterality(*RIGHT),
        **cc(SCT, '266005', 'Structure of lower lobe of right lung'),
    }],
    'subject': PAT,
    'onsetDateTime': '2023-05-25',
    'recordedDate': '2023-05-25',
    'asserter': ONC,
    'note': [{'text': 'EGFR exon 21 L858R-mutated (initial molecular report 2023-06-09).'}],
    'stage': [{
        'summary': cc(SCT, '1222841004', 'IVB (AJCC)', 'Stage IVB (AJCC)'),
        'assessment': [STAGE_GROUP],
        'type': cc(SCT, '897275008', 'American Joint Committee on Cancer, Cancer Staging Manual, 8th edition neoplasm staging system'),
    }],
    'evidence': [{'detail': [
        {'reference': 'DiagnosticReport/mtb-l2345-path-rll-fna'},
        {'reference': 'DiagnosticReport/mtb-l2345-path-left-iliac-biopsy'},
    ]}],
}, [MCODE + 'mcode-primary-cancer-condition'])

related_primary = {'url': 'http://hl7.org/fhir/StructureDefinition/condition-related', 'valueReference': PRIMARY}

add({
    'resourceType': 'Condition', 'id': 'mtb-l2345-bone-metastases',
    'extension': [related_primary],
    **condition_base(),
    'code': {'coding': [
        {'system': SCT, 'code': '94222008', 'display': 'Metastatic malignant neoplasm to bone'},
        {'system': ICD10, 'code': 'C79.51', 'display': 'Secondary malignant neoplasm of bone'},
    ]},
    'bodySite': [
        {'extension': laterality(*LEFT), **cc(SCT, '22356005', 'Bone structure of ilium')},
        cc(SCT, '7610001', 'Bone structure of T10'),
        cc(SCT, '1354727004', 'Bone structure of right sixth rib'),
    ],
    'note': [{'text': 'Left iliac lesion biopsy-proven (2023-05-30).'}],
    'subject': PAT,
    'onsetDateTime': '2023-05-19',
    'recordedDate': '2023-05-30',
    'asserter': ONC,
    'evidence': [{'detail': [
        {'reference': 'DiagnosticReport/mtb-l2345-pet-ct-2023-05-19'},
        {'reference': 'DiagnosticReport/mtb-l2345-path-left-iliac-biopsy'},
    ]}],
}, [MCODE + 'mcode-secondary-cancer-condition'])

add({
    'resourceType': 'Condition', 'id': 'mtb-l2345-lung-metastases',
    'extension': [related_primary],
    **condition_base(),
    'code': cc(SCT, '94391008', 'Metastatic malignant neoplasm to lung'),
    'note': [{'text': 'Bilateral pulmonary nodules, including the contralateral (left) lung.'}],
    'bodySite': [cc(SCT, '44029006', 'Left lung structure'), cc(SCT, '3341006', 'Right lung structure')],
    'subject': PAT,
    'onsetDateTime': '2023-05-15',
    'recordedDate': '2023-05-30',
    'asserter': ONC,
    'evidence': [{'detail': [{'reference': 'DiagnosticReport/mtb-l2345-ct-cap-2023-05-15'}]}],
}, [MCODE + 'mcode-secondary-cancer-condition'])

add({
    'resourceType': 'Condition', 'id': 'mtb-l2345-liver-metastases',
    'extension': [related_primary],
    **condition_base(),
    'code': {'coding': [
        {'system': SCT, 'code': '94381002', 'display': 'Metastatic malignant neoplasm to liver'},
    ]},
    'note': [{'text': 'New lesions in segments 6 and 7 at progression.'}],
    'bodySite': [cc(SCT, '10200004', 'Liver structure')],
    'subject': PAT,
    'onsetDateTime': '2026-09-28',
    'recordedDate': '2026-10-08',
    'asserter': ONC,
    'evidence': [{'detail': [
        {'reference': 'DiagnosticReport/mtb-l2345-ct-cap-2026-09-28'},
        {'reference': 'DiagnosticReport/mtb-l2345-path-liver-biopsy'},
    ]}],
}, [MCODE + 'mcode-secondary-cancer-condition'])

# ---------------------------------------------------------------- Staging (clinical, AJCC 8th)

AJCC8 = cc(SCT, '897275008', 'American Joint Committee on Cancer, Cancer Staging Manual, 8th edition neoplasm staging system')
staging_common = {
    'status': 'final', 'category': obs_cat('imaging', 'Imaging'), 'subject': PAT,
    'focus': [PRIMARY], 'effectiveDateTime': '2023-05-30', 'performer': [ONC], 'method': AJCC8,
}
t = add({'resourceType': 'Observation', 'id': 'mtb-l2345-tnm-clinical-t', **staging_common,
         'code': cc(SCT, '399504009', 'cT category'),
         'valueCodeableConcept': cc(SCT, '1228934000', 'cT2b (AJCC)'), 'note': [{'text': '4.3 cm RLL mass'}]},
        [MCODE + 'mcode-tnm-primary-tumor-category'])
n = add({'resourceType': 'Observation', 'id': 'mtb-l2345-tnm-clinical-n', **staging_common,
         'code': cc(SCT, '399534004', 'cN category'),
         'valueCodeableConcept': cc(SCT, '1229973008', 'cN1 (AJCC)'), 'note': [{'text': 'Ipsilateral (right) hilar adenopathy'}]},
        [MCODE + 'mcode-tnm-regional-nodes-category'])
m = add({'resourceType': 'Observation', 'id': 'mtb-l2345-tnm-clinical-m', **staging_common,
         'code': cc(SCT, '399387003', 'cM category'),
         'valueCodeableConcept': cc(SCT, '1229910003', 'cM1c (AJCC)'),
         'note': [{'text': 'Contralateral pulmonary nodules and multiple extrathoracic (osseous) metastases'}]},
        [MCODE + 'mcode-tnm-distant-metastases-category'])
add({'resourceType': 'Observation', 'id': 'mtb-l2345-tnm-clinical-stage-group', **staging_common,
     'code': cc(SCT, '399537006', 'Clinical TNM stage grouping'),
     'valueCodeableConcept': cc(SCT, '1222841004', 'IVB (AJCC)'),
     'hasMember': [t, n, m]},
    [MCODE + 'mcode-tnm-stage-group'])

# ---------------------------------------------------------------- Tumor size (RLL primary, CT)

TUMOR = add({
    'resourceType': 'BodyStructure', 'id': 'mtb-l2345-tumor-rll',
    'identifier': [{'use': 'usual', 'type': cc('http://hl7.org/fhir/resource-types', 'BodyStructure', 'BodyStructure'),
                    'system': 'http://elimu.io/fhir/tumor-id', 'value': 'L2345-RLL-primary'}],
    'morphology': cc(SCT, '1187332001', 'Adenocarcinoma'),
    'location': cc(SCT, '266005', 'Structure of lower lobe of right lung'),
    'locationQualifier': [cc(SCT, '24028007', 'Right')],
    'description': 'Primary tumor, right lower lobe of lung',
    'patient': PAT,
}, [MCODE + 'mcode-tumor'])


def tumor_size(rid, date, longest, other):
    add({'resourceType': 'Observation', 'id': rid, 'status': 'final',
         'category': obs_cat('imaging', 'Imaging'),
         'code': cc(LOINC, '21889-1', 'Size Tumor'),
         'subject': PAT, 'focus': [TUMOR], 'effectiveDateTime': date, 'performer': [RAD],
         'method': cc(SCT, '77477000', 'Computerized axial tomography'),
         'component': [
             {'code': cc(LOINC, '33728-7', 'Size.maximum dimension in Tumor'), 'valueQuantity': qty(longest, 'cm')},
             {'code': cc(LOINC, '33729-5', 'Size additional dimension in Tumor'), 'valueQuantity': qty(other, 'cm')},
         ]}, [MCODE + 'mcode-tumor-size'])


tumor_size('mtb-l2345-tumor-size-2023-05-15', '2023-05-15', 4.3, 3.6)
tumor_size('mtb-l2345-tumor-size-2024-03-11', '2024-03-11', 2.1, 1.7)
tumor_size('mtb-l2345-tumor-size-2026-09-28', '2026-09-28', 3.4, 2.8)

# ---------------------------------------------------------------- Imaging & pathology reports

DR_CAT_RAD = cc(LOINC, 'LP29684-5', 'Radiology')
DR_CAT_PATH = cc(LOINC, 'LP7839-6', 'Pathology')
DR_CAT_GE = cc('http://terminology.hl7.org/CodeSystem/v2-0074', 'GE', 'Genetics')


def report(rid, date, cat, code, conclusion, performer, specimen=None, conclusion_code=None, results=None,
           genomics=False):
    # A genomic profiling report with structured variants is also an mCODE genomics report
    r = {'resourceType': 'DiagnosticReport', 'id': rid, 'status': 'final', 'category': [cat], 'code': code,
         'subject': PAT, 'effectiveDateTime': date, 'issued': f'{date}T17:00:00-04:00',
         'performer': [performer], 'conclusion': conclusion}
    profiles = [USCORE + 'us-core-diagnosticreport-note']
    if genomics:
        r['category'].append(DR_CAT_GE)
        profiles.insert(0, MCODE + 'mcode-genomics-report')
    if specimen:
        r['specimen'] = [specimen]
    if conclusion_code:
        r['conclusionCode'] = [conclusion_code]
    if results:
        r['result'] = results
    add(r, profiles)


# Genomic variants and therapeutic implications (incorporating-curated-lab-report-design.md §2.4): variants in
# the reportable genes from the genomics operations server's L2345 data (GRCh38; recorded as somatic, no VAF),
# and therapeutic implications in the server's CIViC shape
GR_TBD = 'http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs'
HGVS = 'http://varnomen.hgvs.org'
LAB_GE = obs_cat('laboratory', 'Laboratory') + [DR_CAT_GE]
NSCLC = cc('https://disease-ontology.org', 'DOID:3908', 'Lung Non-small Cell Carcinoma')
LUNG_CANCER = cc('https://disease-ontology.org', 'DOID:1324', 'Lung Cancer')
LUNG_ADENO = cc('https://disease-ontology.org', 'DOID:3910', 'Lung Adenocarcinoma')
LEVEL_A = 'A (Validated association)'
LEVEL_B = 'B (Clinical evidence)'


def component(code, display, value, system=LOINC):
    value = value if any(key.startswith('value') for key in value) else {'valueCodeableConcept': value}
    return {'code': cc(system, code, display), **value}


def genomic_variant(rid, date, specimen, gene, hgnc, c_hgvs, p_hgvs, consequence, refseq, spdi, ref, alt, start):
    return add({
        'resourceType': 'Observation', 'id': rid, 'status': 'final', 'category': LAB_GE,
        'code': cc(LOINC, '69548-6', 'Genetic variant assessment'), 'subject': PAT,
        'effectiveDateTime': date, 'issued': f'{date}T17:00:00-04:00', 'performer': [PATH],
        'valueCodeableConcept': cc(LOINC, 'LA9633-4', 'Present'),
        'method': cc(LOINC, 'LA26398-0', 'Sequencing'), 'specimen': specimen,
        'component': [
            component('48018-6', 'Gene studied [ID]', cc('http://www.genenames.org', hgnc, gene)),
            component('48004-6', 'DNA change (c.HGVS)', {'coding': [{'system': HGVS, 'code': c_hgvs}]}),
            component('48005-3', 'Amino acid change (pHGVS)', {'coding': [{'system': HGVS, 'code': p_hgvs}]}),
            component('molecular-consequence', 'Molecular Consequence',
                      cc('http://www.sequenceontology.org/', *consequence), system=GR_TBD),
            component('48002-0', 'Genomic source class [Type]', cc(LOINC, 'LA6684-0', 'Somatic')),
            component('48013-7', 'Genomic reference sequence [ID]',
                      {'coding': [{'system': 'http://www.ncbi.nlm.nih.gov/nuccore', 'code': refseq}]}),
            component('81252-9', 'Discrete genetic variant',
                      cc('https://api.ncbi.nlm.nih.gov/variation/v0/', spdi, spdi)),
            component('69547-8', 'Genomic ref allele [ID]', {'valueString': ref}),
            # A pure deletion has no alt allele (FHIR strings can't be empty); the server leaves it out too
            *([component('69551-0', 'Genomic alt allele [ID]', {'valueString': alt})] if alt else []),
            component('92822-6', 'Genomic coordinate system [Type]', cc(LOINC, 'LA30100-4', '0-based interval counting')),
            component('81254-5', 'Genomic allele start-end', {'valueRange': {'low': {'value': start}}}),
        ],
    }, [MCODE + 'mcode-genomic-variant'])


def tx_implication(rid, date, variant, rxcui, drug, implication, evidence, phenotype):
    return add({
        'resourceType': 'Observation', 'id': rid, 'status': 'final', 'category': obs_cat('laboratory', 'Laboratory'),
        'code': cc(GR_TBD, 'therapeutic-implication', 'Therapeutic Implication'), 'subject': PAT,
        'effectiveDateTime': date, 'issued': f'{date}T17:00:00-04:00', 'performer': [PATH],
        'derivedFrom': [variant],
        'component': [
            component('81259-4', 'Associated phenotype', phenotype),
            component('93044-6', 'Level of evidence', {'text': evidence}),
            component('51963-7', 'Medication assessed [ID]', cc(RXN, rxcui, drug)),
            component('predicted-therapeutic-implication', 'Predicted Therapeutic Implication', {'text': implication},
                      system=GR_TBD),
        ],
    }, ['http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication'])


MISSENSE = ('SO:0001583', 'missense_variant')
FRAMESHIFT = ('SO:0001589', 'frameshift_variant')
STOP_GAINED = ('SO:0001587', 'stop_gained')


def egfr_l858r(year, date, specimen):
    # EGFR L858R is in both profiling reports; each report's copy carries its own implications
    variant = genomic_variant(f'mtb-l2345-genomic-variant-egfr-l858r-{year}', date, specimen, 'EGFR', 'HGNC:3236',
                              'NM_005228.5:c.2573T>G', 'NP_005219.2:p.Leu858Arg', MISSENSE, 'NC_000007.14',
                              'NC_000007.14:55191821:T:G', 'T', 'G', 55191821)
    return [variant,
            # CIViC EID2994, EID2997
            tx_implication(f'mtb-l2345-tx-implication-erlotinib-l858r-{year}', date, variant, '337525', 'Erlotinib',
                           'Supports Sensitivity/Response', LEVEL_A, NSCLC),
            tx_implication(f'mtb-l2345-tx-implication-afatinib-l858r-{year}', date, variant, '1430438', 'Afatinib',
                           'Supports Sensitivity/Response', LEVEL_A, NSCLC)]


report('mtb-l2345-cxr-2023-05-10', '2023-05-10', DR_CAT_RAD, cc(LOINC, '36643-5', 'XR Chest 2 Views'),
       '4 cm right lower lobe mass. Right hilar fullness. No pleural effusion. Recommend CT chest.', RAD)
report('mtb-l2345-ct-cap-2023-05-15', '2023-05-15', DR_CAT_RAD,
       cc(LOINC, '72254-6', 'CT Chest and Abdomen and Pelvis W contrast IV'),
       'Spiculated 4.3 x 3.6 cm right lower lobe mass. Enlarged right hilar lymph nodes (short axis up to 1.6 cm); '
       'no mediastinal or contralateral hilar adenopathy. Multiple bilateral pulmonary nodules up to 0.9 cm, including '
       'four in the left lung, suspicious for metastases. Lytic lesions in the left ilium (3.1 cm), T10 vertebral body '
       'and right 6th rib. Liver, adrenals unremarkable. No pleural or pericardial effusion.', RAD,
       results=[{'reference': 'Observation/mtb-l2345-tumor-size-2023-05-15'}])
report('mtb-l2345-pet-ct-2023-05-19', '2023-05-19', DR_CAT_RAD,
       cc(LOINC, '81554-8', 'PET+CT Guidance for localization of tumor of Skull base to mid-thigh-- W 18F-FDG IV'),
       'FDG-avid right lower lobe mass (SUVmax 12.4) and right hilar nodes (SUVmax 6.8). Bilateral FDG-avid pulmonary '
       'nodules. FDG-avid osseous metastases: left ilium (SUVmax 9.7), T10, right 6th rib. No FDG-avid mediastinal, '
       'hepatic or adrenal disease.', RAD)
report('mtb-l2345-mri-brain-2023-05-22', '2023-05-22', DR_CAT_RAD, cc(LOINC, '24587-8', 'MR Brain WO and W contrast IV'),
       'No evidence of intracranial metastatic disease.', RAD)
report('mtb-l2345-path-rll-fna', '2023-05-25', DR_CAT_PATH, cc(LOINC, '11526-1', 'Pathology study'),
       'Right lower lobe mass, CT-guided FNA: non-small cell carcinoma, consistent with adenocarcinoma. '
       'Immunostains on cell block: TTF-1 positive, Napsin A positive, p40 negative. Cell block adequate for '
       'ancillary studies (PD-L1 and molecular testing sent).', PATH, SPEC_LUNG,
       cc(SCT, '1187332001', 'Adenocarcinoma'))
report('mtb-l2345-path-left-iliac-biopsy', '2023-05-30', DR_CAT_PATH, cc(LOINC, '11526-1', 'Pathology study'),
       'Left iliac bone lesion, CT-guided core biopsy: metastatic adenocarcinoma, TTF-1 positive, morphologically '
       'similar to the right lower lobe FNA, consistent with a lung primary.', PATH, SPEC_ILIUM,
       cc(SCT, '1187332001', 'Adenocarcinoma'))
report('mtb-l2345-molecular-2023-06-09', '2023-06-09', DR_CAT_PATH, cc(LOINC, '51969-4', 'Genetic analysis report'),
       'Initial tumor genomic profiling (RLL FNA cell block, targeted NGS): EGFR exon 21 L858R detected. '
       'No ALK, ROS1, RET or NTRK fusions; no MET exon 14 skipping; KRAS and BRAF wild type. EGFR L858R: predicted '
       'sensitivity to erlotinib and afatinib in non-small cell lung cancer (level A).', PATH, SPEC_LUNG,
       results=egfr_l858r(2023, '2023-06-09', SPEC_LUNG), genomics=True)
report('mtb-l2345-ct-cap-2024-03-11', '2024-03-11', DR_CAT_RAD,
       cc(LOINC, '72254-6', 'CT Chest and Abdomen and Pelvis W contrast IV'),
       'Continued decrease in right lower lobe mass, now 2.1 x 1.7 cm (nadir). Right hilar nodes normal in size. '
       'Pulmonary nodules resolved or decreased. Sclerotic (treated) osseous metastases. No new disease.', RAD,
       results=[{'reference': 'Observation/mtb-l2345-tumor-size-2024-03-11'}])
report('mtb-l2345-ct-cap-2026-09-28', '2026-09-28', DR_CAT_RAD,
       cc(LOINC, '72254-6', 'CT Chest and Abdomen and Pelvis W contrast IV'),
       'Progression of disease. Right lower lobe mass enlarged to 3.4 x 2.8 cm (from nadir 2.1 x 1.7 cm). '
       'Two new hypoattenuating liver lesions in segments 6 (2.4 cm) and 7 (1.5 cm), consistent with metastases. '
       'New 1.1 cm left lower lobe nodule. Osseous metastases stable and sclerotic. No new bone lesions.', RAD,
       results=[{'reference': 'Observation/mtb-l2345-tumor-size-2026-09-28'}])
report('mtb-l2345-mri-brain-2026-10-02', '2026-10-02', DR_CAT_RAD, cc(LOINC, '24587-8', 'MR Brain WO and W contrast IV'),
       'No evidence of intracranial metastatic disease.', RAD)
report('mtb-l2345-path-liver-biopsy', '2026-10-06', DR_CAT_PATH, cc(LOINC, '11526-1', 'Pathology study'),
       'Liver segment 6 lesion, CT-guided core biopsy: metastatic adenocarcinoma, TTF-1 positive, consistent with '
       'the known lung primary. No small cell or squamous transformation identified. Tissue sent for repeat '
       'comprehensive genomic profiling.', PATH, SPEC_LIVER,
       cc(SCT, '1187332001', 'Adenocarcinoma'))

# Repeat profiling at progression: the profile the genomics server holds
REPEAT_DATE = '2026-10-16'


def repeat_variant(slug, gene, hgnc, c_hgvs, p_hgvs, consequence, refseq, spdi, ref, alt, start):
    return genomic_variant(f'mtb-l2345-genomic-variant-{slug}', REPEAT_DATE, SPEC_LIVER, gene, hgnc, c_hgvs, p_hgvs,
                           consequence, refseq, spdi, ref, alt, start)


def repeat_implication(slug, variant, rxcui, drug, implication, evidence, phenotype=NSCLC):
    return tx_implication(f'mtb-l2345-tx-implication-{slug}', REPEAT_DATE, variant, rxcui, drug, implication,
                          evidence, phenotype)


EGFR_T790M = repeat_variant('egfr-t790m', 'EGFR', 'HGNC:3236', 'NM_005228.5:c.2369C>T', 'NP_005219.2:p.Thr790Met',
                            MISSENSE, 'NC_000007.14', 'NC_000007.14:55181377:C:T', 'C', 'T', 55181377)
# The server gives only p.Lys745fs; the full form is from Mutalyzer
EGFR_K745FS = repeat_variant('egfr-k745fs', 'EGFR', 'HGNC:3236', 'NM_005228.5:c.2232_2250del',
                             'NP_005219.2:p.(Lys745HisfsTer15)', FRAMESHIFT, 'NC_000007.14',
                             'NC_000007.14:55174768:CAAGGAATTAAGAGAAGCA:', 'CAAGGAATTAAGAGAAGCA', '', 55174768)
KRAS_G12C = repeat_variant('kras-g12c', 'KRAS', 'HGNC:6407', 'NM_004985.5:c.34G>T', 'NP_004976.2:p.Gly12Cys',
                           MISSENSE, 'NC_000012.12', 'NC_000012.12:25245350:C:A', 'C', 'A', 25245350)
KRAS_G12A = repeat_variant('kras-g12a', 'KRAS', 'HGNC:6407', 'NM_004985.5:c.35G>C', 'NP_004976.2:p.Gly12Ala',
                           MISSENSE, 'NC_000012.12', 'NC_000012.12:25245349:C:G', 'C', 'G', 25245349)
BRAF_V600E = repeat_variant('braf-v600e', 'BRAF', 'HGNC:1097', 'NM_004333.6:c.1799T>A', 'NP_004324.2:p.Val600Glu',
                            MISSENSE, 'NC_000007.14', 'NC_000007.14:140753335:A:T', 'A', 'T', 140753335)
ALK_F1245V = repeat_variant('alk-f1245v', 'ALK', 'HGNC:427', 'NM_004304.5:c.3733T>G', 'NP_004295.2:p.Phe1245Val',
                            MISSENSE, 'NC_000002.12', 'NC_000002.12:29213993:A:C', 'A', 'C', 29213993)
BRCA2_Q66X = repeat_variant('brca2-q66x', 'BRCA2', 'HGNC:1101', 'NM_000059.4:c.196C>T', 'NP_000050.3:p.Gln66Ter',
                            STOP_GAINED, 'NC_000013.11', 'NC_000013.11:32319204:C:T', 'C', 'T', 32319204)
# The server gives no c./p. change; both are from Mutalyzer
STK11_FS = repeat_variant('stk11-k48fs', 'STK11', 'HGNC:11389', 'NM_000455.5:c.143_158del',
                          'NP_000446.1:p.(Lys48ThrfsTer11)', FRAMESHIFT, 'NC_000019.10',
                          'NC_000019.10:1207054:AAGTACCTGATGGGGGA:A', 'AAGTACCTGATGGGGGA', 'A', 1207054)

REPEAT_RESULTS = egfr_l858r(2026, REPEAT_DATE, SPEC_LIVER) + [
    EGFR_T790M,
    # CIViC EID1592, EID238
    repeat_implication('osimertinib-t790m', EGFR_T790M, '1721560', 'Osimertinib', 'Supports Sensitivity/Response', LEVEL_A),
    repeat_implication('erlotinib-t790m', EGFR_T790M, '337525', 'Erlotinib', 'Supports Resistance', LEVEL_A),
    EGFR_K745FS,
    KRAS_G12C,
    # CIViC EID9431 (not in the genomics server's CIViC subset), EID2257
    repeat_implication('sotorasib-g12c', KRAS_G12C, '2550714', 'Sotorasib', 'Supports Sensitivity/Response', LEVEL_A),
    repeat_implication('gefitinib-g12c', KRAS_G12C, '328134', 'Gefitinib', 'Supports Resistance', LEVEL_B, LUNG_CANCER),
    KRAS_G12A,
    # CIViC EID251
    repeat_implication('erlotinib-g12a', KRAS_G12A, '337525', 'Erlotinib', 'Supports Resistance', LEVEL_B, LUNG_ADENO),
    BRAF_V600E,
    # CIViC EID3017 (in practice dabrafenib + trametinib), EID5958
    repeat_implication('trametinib-v600e', BRAF_V600E, '1425099', 'Trametinib', 'Supports Sensitivity/Response', LEVEL_A),
    repeat_implication('vemurafenib-v600e', BRAF_V600E, '1147220', 'Vemurafenib', 'Supports Sensitivity/Response', LEVEL_B),
    ALK_F1245V, BRCA2_Q66X, STK11_FS,
]

report('mtb-l2345-molecular-2026-10-16', REPEAT_DATE, DR_CAT_PATH, cc(LOINC, '51969-4', 'Genetic analysis report'),
       'Repeat comprehensive genomic profiling at progression (liver segment 6 core biopsy, targeted NGS): EGFR L858R '
       '(known) with EGFR T790M and EGFR c.2232_2250del (p.Lys745HisfsTer15); KRAS G12C and G12A; BRAF V600E; '
       'ALK F1245V (point mutation; no ALK fusion); BRCA2 Q66*; STK11 c.143_158del (p.Lys48ThrfsTer11). '
       'Predicted, in non-small cell lung cancer: EGFR T790M, sensitivity to osimertinib and resistance to erlotinib '
       '(level A); KRAS G12C, sensitivity to sotorasib (level A); BRAF V600E, sensitivity to trametinib (level A) and '
       'vemurafenib (level B); KRAS G12C and G12A, resistance to gefitinib and erlotinib (level B).',
       PATH, SPEC_LIVER, results=REPEAT_RESULTS, genomics=True)

# ---------------------------------------------------------------- Tumor markers


def tumor_marker(rid, date, code, value, specimen=None, interpretation=None, note=None):
    r = {'resourceType': 'Observation', 'id': rid, 'status': 'final',
         'extension': [ext('mcode-related-condition', 'valueReference', PRIMARY)],
         'category': obs_cat('laboratory', 'Laboratory'), 'code': code, 'subject': PAT,
         'effectiveDateTime': date, 'performer': [PATH], **value}
    if specimen:
        r['specimen'] = specimen
    if interpretation:
        r['interpretation'] = [interpretation]
    if note:
        r['note'] = [{'text': note}]
    add(r, [MCODE + 'mcode-tumor-marker-test'])


tumor_marker('mtb-l2345-pdl1-tps-2023-06-05', '2023-06-05',
             cc(LOINC, '83053-9', 'Cells.programmed cell death ligand 1/Viable tumor cells in Tissue by Immune stain',
                'PD-L1 tumor proportion score'),
             {'valueQuantity': qty(5, '%')}, SPEC_LUNG,
             note='PD-L1 IHC clone 22C3 pharmDx on RLL FNA cell block. TPS 5% (low expression, 1-49%).')

HIGH = cc('http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation', 'H', 'High')
for date, value in [('2023-05-25', 48.2), ('2023-09-11', 6.4), ('2024-03-11', 3.1), ('2025-03-10', 3.4),
                    ('2026-03-16', 4.0), ('2026-09-28', 21.7)]:
    tumor_marker(f'mtb-l2345-cea-{date}', date,
                 cc(LOINC, '2039-6', 'Carcinoembryonic Ag [Mass/volume] in Serum or Plasma'),
                 {'valueQuantity': qty(value, 'ng/mL'),
                  'referenceRange': [{'high': qty(3.0, 'ng/mL'), 'text': '<3.0 ng/mL (non-smoker)'}]},
                 interpretation=HIGH if value > 3.0 else None)

# ---------------------------------------------------------------- Disease status

for date, code, display, note in [
    ('2023-09-11', '268910001', "Patient's condition improved", 'Partial response after 4 cycles of carboplatin/pemetrexed + osimertinib.'),
    ('2024-03-11', '268910001', "Patient's condition improved", 'Continued response; RLL mass at nadir.'),
    ('2025-03-10', '359746009', "Patient's condition stable", 'Stable disease on osimertinib + pemetrexed maintenance.'),
    ('2026-03-16', '359746009', "Patient's condition stable", 'Stable disease on osimertinib + pemetrexed maintenance.'),
    ('2026-09-28', '271299001', "Patient's condition worsened", 'Progression: enlarging RLL primary, new liver metastases, new LLL nodule.'),
]:
    add({'resourceType': 'Observation', 'id': f'mtb-l2345-disease-status-{date}',
         'extension': [ext('mcode-cancer-disease-status-evidence-type', 'valueCodeableConcept',
                           cc(SCT, '363679005', 'Imaging'))],
         'status': 'final', 'category': obs_cat('exam', 'Exam'),
         'code': cc(LOINC, '97509-4', 'Cancer disease progression'),
         'subject': PAT, 'focus': [PRIMARY], 'effectiveDateTime': date, 'performer': [ONC],
         'valueCodeableConcept': cc(SCT, code, display), 'note': [{'text': note}]},
        [MCODE + 'mcode-cancer-disease-status'])

# ---------------------------------------------------------------- ECOG

for date, score, la, display in [
    ('2023-05-25', 1, 'LA9623-5', 'Restricted in physically strenuous activity but ambulatory and able to carry out work of a light or sedentary nature, e.g., light house work, office work'),
    ('2023-09-11', 0, 'LA9622-7', 'Fully active, able to carry on all pre-disease performance without restriction'),
    ('2026-10-08', 1, 'LA9623-5', 'Restricted in physically strenuous activity but ambulatory and able to carry out work of a light or sedentary nature, e.g., light house work, office work'),
]:
    add({'resourceType': 'Observation', 'id': f'mtb-l2345-ecog-{date}', 'status': 'final',
         'category': obs_cat('survey', 'Survey'),
         'code': cc(LOINC, '89247-1', 'ECOG Performance Status score'),
         'subject': PAT, 'effectiveDateTime': date, 'performer': [ONC], 'valueInteger': score,
         'interpretation': [cc(LOINC, la, display)]},
        [MCODE + 'mcode-ecog-performance-status'])

# ---------------------------------------------------------------- Cancer-related medications

PALLIATIVE = ext('mcode-procedure-intent', 'valueCodeableConcept', cc(SCT, '363676003', 'Palliative intent'))
SUPPORTIVE = ext('mcode-procedure-intent', 'valueCodeableConcept', cc(SCT, '399707004', 'Supportive - procedure intent'))
NO_RESPONSE = cc(SCT, '266721009', 'Absent response to treatment')


def med_request(rid, status, rxcui, display, start, dosage_text, end=None, status_reason=None, intent_ext=PALLIATIVE,
                note=None, reason=PRIMARY, timing=None, dose=None, route=None):
    bounds = {'start': start, **({'end': end} if end else {})}
    dosage = {'text': dosage_text, 'timing': {'repeat': {'boundsPeriod': bounds, **(timing or {})}}}
    if route:
        dosage['route'] = route
    if dose:
        dosage['doseAndRate'] = [dose]
    r = {'resourceType': 'MedicationRequest', 'id': rid, 'extension': [intent_ext], 'status': status,
         'intent': 'order', 'category': [cc('http://terminology.hl7.org/CodeSystem/medicationrequest-category',
                                            'outpatient', 'Outpatient')],
         'medicationCodeableConcept': cc(RXN, rxcui, display), 'subject': PAT, 'authoredOn': start,
         'requester': ONC, 'reasonReference': [reason], 'dosageInstruction': [dosage]}
    if status_reason:
        r['statusReason'] = status_reason
    if note:
        r['note'] = [{'text': note}]
    add(r, [MCODE + 'mcode-cancer-related-medication-request'])


ORAL = cc(SCT, '26643006', 'Oral route')
IV = cc(SCT, '47625008', 'Intravenous route')

med_request('mtb-l2345-medrx-osimertinib', 'active', '1721581', 'osimertinib 80 MG Oral Tablet', '2023-06-19',
            'Osimertinib 80 mg by mouth once daily', timing={'frequency': 1, 'period': 1, 'periodUnit': 'd'},
            dose={'doseQuantity': qty(80, 'mg')}, route=ORAL,
            note='First-line, with carboplatin/pemetrexed (FLAURA2 regimen) for EGFR L858R-mutated stage IVB adenocarcinoma. '
                 'Continued after the 2026-09-28 progression scan pending repeat genomic profiling of the liver lesion '
                 'and molecular tumor board review.')
med_request('mtb-l2345-medrx-carboplatin', 'completed', '597195', 'carboplatin 10 MG/ML Injectable Solution',
            '2023-06-19', 'Carboplatin AUC 5 IV on day 1 every 21 days x 4 cycles', end='2023-08-21',
            timing={'count': 4, 'frequency': 1, 'period': 21, 'periodUnit': 'd'}, route=IV,
            note='Cycles given 2023-06-19, 2023-07-10, 2023-07-31, 2023-08-21.')
med_request('mtb-l2345-medrx-pemetrexed', 'stopped', '1728072', 'pemetrexed 500 MG Injection', '2023-06-19',
            'Pemetrexed 500 mg/m2 IV on day 1 every 21 days (4 cycles with carboplatin, then maintenance)',
            end='2026-09-28', status_reason=NO_RESPONSE,
            timing={'frequency': 1, 'period': 21, 'periodUnit': 'd'}, route=IV,
            dose={'doseQuantity': qty(500, 'mg/m2', 'mg/m2')},
            note='With folic acid and vitamin B12 supplementation. Maintenance stopped at progression on 2026-09-28.')
med_request('mtb-l2345-medrx-zoledronic-acid', 'active', '1865665', 'zoledronic acid 4 MG Injection', '2023-07-10',
            'Zoledronic acid 4 mg IV every 12 weeks', intent_ext=SUPPORTIVE,
            timing={'frequency': 1, 'period': 12, 'periodUnit': 'wk'}, route=IV,
            reason={'reference': 'Condition/mtb-l2345-bone-metastases', 'display': 'Bone metastases'},
            note='Bone-modifying agent for osseous metastases. Calcium and vitamin D supplementation.')

# ---------------------------------------------------------------- Radiotherapy (palliative, left ilium)

RT_VOLUME = add({
    'resourceType': 'BodyStructure', 'id': 'mtb-l2345-rt-volume-left-ilium',
    'identifier': [{'use': 'usual', 'system': 'http://elimu.io/fhir/rt-volume', 'value': 'PTV_L_Ilium'}],
    'morphology': cc(SCT, '228793007', 'Planning target volume'),
    'location': cc(SCT, '22356005', 'Bone structure of ilium'),
    'locationQualifier': [cc(SCT, '7771000', 'Left')],
    'description': 'Left iliac osseous metastasis',
    'patient': PAT,
}, [MCODE + 'mcode-radiotherapy-volume'])

add({
    'resourceType': 'Procedure', 'id': 'mtb-l2345-rt-course-left-ilium',
    'extension': [
        PALLIATIVE,
        {'url': MCODE + 'mcode-radiotherapy-modality-and-technique', 'extension': [
            {'url': MCODE + 'mcode-radiotherapy-modality', 'valueCodeableConcept':
                cc(SCT, '1156506007', 'External beam radiation therapy using photons')},
            {'url': MCODE + 'mcode-radiotherapy-technique', 'valueCodeableConcept':
                cc(SCT, '1162782007', 'Three dimensional external beam radiation therapy')},
        ]},
        {'url': MCODE + 'mcode-radiotherapy-sessions', 'valueUnsignedInt': 5},
        {'url': MCODE + 'mcode-radiotherapy-dose-delivered-to-volume', 'extension': [
            {'url': 'volume', 'valueReference': RT_VOLUME},
            {'url': 'totalDoseDelivered', 'valueQuantity': qty(2000, 'cGy')},
            {'url': 'fractionsDelivered', 'valueUnsignedInt': 5},
        ]},
    ],
    'status': 'completed',
    'code': cc(SCT, '1217123003', 'Radiotherapy course of treatment'),
    'subject': PAT,
    'performedPeriod': {'start': '2023-06-05', 'end': '2023-06-09'},
    'performer': [{'actor': RADONC}],
    'reasonReference': [{'reference': 'Condition/mtb-l2345-bone-metastases', 'display': 'Bone metastases'}],
    'bodySite': [{'extension': laterality(*LEFT), **cc(SCT, '22356005', 'Bone structure of ilium')}],
    'note': [{'text': 'Palliative RT for painful left iliac metastasis, 20 Gy in 5 fractions. Pain resolved.'}],
}, [MCODE + 'mcode-radiotherapy-course-summary'])

# ---------------------------------------------------------------- Diagnostic procedures (biopsies)


def biopsy(rid, date, code, site, report_id, reason=PRIMARY):
    add({'resourceType': 'Procedure', 'id': rid,
         'status': 'completed', 'code': code, 'subject': PAT, 'performedDateTime': date,
         'performer': [{'actor': RAD}], 'reasonReference': [reason], 'bodySite': [site],
         'report': [{'reference': f'DiagnosticReport/{report_id}'}]},
        [USCORE + 'us-core-procedure'])


biopsy('mtb-l2345-proc-rll-fna', '2023-05-25', cc(SCT, '432231006', 'Fine needle aspiration biopsy of lung'),
       {'extension': laterality(*RIGHT), **cc(SCT, '266005', 'Structure of lower lobe of right lung')},
       'mtb-l2345-path-rll-fna')
biopsy('mtb-l2345-proc-left-iliac-biopsy', '2023-05-30', cc(SCT, '432371004', 'CT guided biopsy of bone'),
       {'extension': laterality(*LEFT), **cc(SCT, '22356005', 'Bone structure of ilium')},
       'mtb-l2345-path-left-iliac-biopsy')
biopsy('mtb-l2345-proc-liver-biopsy', '2026-10-06',
       cc(SCT, '407971000119109', 'CT guided percutaneous needle biopsy of liver'),
       cc(SCT, '10200004', 'Liver structure'), 'mtb-l2345-path-liver-biopsy')

# ---------------------------------------------------------------- Comorbidities, allergies, other meds

htn = add({'resourceType': 'Condition', 'id': 'mtb-l2345-hypertension', **condition_base(),
           'code': cc(SCT, '59621000', 'Essential hypertension'),
           'subject': PAT, 'onsetDateTime': '2014', 'recordedDate': '2014-08-20'},
          [USCORE + 'us-core-condition-problems-health-concerns'])
hld = add({'resourceType': 'Condition', 'id': 'mtb-l2345-hyperlipidemia', **condition_base(),
           'code': cc(SCT, '55822004', 'Hyperlipidemia'),
           'subject': PAT, 'onsetDateTime': '2011', 'recordedDate': '2011-02-10'},
          [USCORE + 'us-core-condition-problems-health-concerns'])

add({'resourceType': 'Observation', 'id': 'mtb-l2345-comorbidities',
     'extension': [ext('mcode-related-condition', 'valueReference', htn),
                   ext('mcode-related-condition', 'valueReference', hld)],
     'status': 'final', 'category': obs_cat('exam', 'Exam'),
     'code': cc(SCT, '398192003', 'Co-morbid conditions'), 'subject': PAT, 'focus': [PRIMARY],
     'effectiveDateTime': '2026-10-08', 'performer': [ONC]},
    [MCODE + 'mcode-comorbidities'])

for rid, rxcui, display, text, reason in [
    ('mtb-l2345-medrx-lisinopril', '314077', 'lisinopril 20 MG Oral Tablet', 'Lisinopril 20 mg by mouth once daily', htn),
    ('mtb-l2345-medrx-atorvastatin', '617311', 'atorvastatin 40 MG Oral Tablet', 'Atorvastatin 40 mg by mouth once daily', hld),
]:
    add({'resourceType': 'MedicationRequest', 'id': rid, 'status': 'active', 'intent': 'order',
         'medicationCodeableConcept': cc(RXN, rxcui, display), 'subject': PAT, 'authoredOn': '2026-03-16',
         'requester': ONC, 'reasonReference': [reason],
         'dosageInstruction': [{'text': text, 'route': ORAL,
                                'timing': {'repeat': {'frequency': 1, 'period': 1, 'periodUnit': 'd'}}}]},
        [USCORE + 'us-core-medicationrequest'])

add({'resourceType': 'AllergyIntolerance', 'id': 'mtb-l2345-allergy-penicillin',
     'clinicalStatus': cc('http://terminology.hl7.org/CodeSystem/allergyintolerance-clinical', 'active', 'Active'),
     'verificationStatus': cc('http://terminology.hl7.org/CodeSystem/allergyintolerance-verification', 'confirmed', 'Confirmed'),
     'type': 'allergy', 'category': ['medication'], 'criticality': 'low',
     'code': cc(RXN, '7980', 'penicillin G', 'Penicillin G'),
     'patient': PAT, 'recordedDate': '2010-01-05',
     'reaction': [{'manifestation': [cc(SCT, '271807003', 'Eruption', 'Rash')], 'severity': 'mild'}]},
    [USCORE + 'us-core-allergyintolerance'])

add({'resourceType': 'Observation', 'id': 'mtb-l2345-smoking-status', 'status': 'final',
     'category': obs_cat('social-history', 'Social History'),
     'code': cc(LOINC, '72166-2', 'Tobacco smoking status'), 'subject': PAT,
     'effectiveDateTime': '2023-05-10', 'performer': [ONC],
     'valueCodeableConcept': cc(SCT, '266919005', 'Never smoked tobacco')},
    [USCORE + 'us-core-smokingstatus'])

# ---------------------------------------------------------------- Vital signs

VS_CAT = obs_cat('vital-signs', 'Vital Signs')
for date, weight in [('2023-05-25', 78.0), ('2026-10-08', 74.2)]:
    add({'resourceType': 'Observation', 'id': f'mtb-l2345-weight-{date}', 'status': 'final', 'category': VS_CAT,
         'code': cc(LOINC, '29463-7', 'Body weight'), 'subject': PAT, 'effectiveDateTime': date, 'performer': [ONC],
         'valueQuantity': qty(weight, 'kg')},
        ['http://hl7.org/fhir/StructureDefinition/bodyweight', USCORE + 'us-core-body-weight'])
add({'resourceType': 'Observation', 'id': 'mtb-l2345-height-2023-05-25', 'status': 'final', 'category': VS_CAT,
     'code': cc(LOINC, '8302-2', 'Body height'), 'subject': PAT, 'effectiveDateTime': '2023-05-25', 'performer': [ONC],
     'valueQuantity': qty(175, 'cm')},
    ['http://hl7.org/fhir/StructureDefinition/bodyheight', USCORE + 'us-core-body-height'])

# ---------------------------------------------------------------- Current labs (2026-10-08)

LOW = cc('http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation', 'L', 'Low')
for code, display, value, unit, low, high, interp in [
    ('2160-0', 'Creatinine [Mass/volume] in Serum or Plasma', 1.1, 'mg/dL', 0.7, 1.3, None),
    ('98979-8', 'Glomerular filtration rate [Volume Rate/Area] in Serum, Plasma or Blood by Creatinine-based formula (CKD-EPI 2021)/1.73 sq M', 72, 'mL/min/{1.73_m2}', 60, None, None),
    ('1920-8', 'Aspartate aminotransferase [Enzymatic activity/volume] in Serum or Plasma', 58, 'U/L', 10, 40, HIGH),
    ('1742-6', 'Alanine aminotransferase [Enzymatic activity/volume] in Serum or Plasma', 63, 'U/L', 7, 56, HIGH),
    ('1975-2', 'Bilirubin.total [Mass/volume] in Serum or Plasma', 0.9, 'mg/dL', 0.1, 1.2, None),
    ('718-7', 'Hemoglobin [Mass/volume] in Blood', 12.1, 'g/dL', 13.5, 17.5, LOW),
    ('6690-2', 'Leukocytes [#/volume] in Blood by Automated count', 5.8, '10*3/uL', 4.0, 11.0, None),
    ('751-8', 'Neutrophils [#/volume] in Blood by Automated count', 3.6, '10*3/uL', 1.8, 7.7, None),
    ('777-3', 'Platelets [#/volume] in Blood by Automated count', 212, '10*3/uL', 150, 400, None),
]:
    rr = {**({'low': qty(low, unit)} if low is not None else {}), **({'high': qty(high, unit)} if high is not None else {})}
    r = {'resourceType': 'Observation', 'id': f'mtb-l2345-lab-{code}-2026-10-08', 'status': 'final',
         'category': obs_cat('laboratory', 'Laboratory'), 'code': cc(LOINC, code, display), 'subject': PAT,
         'effectiveDateTime': '2026-10-08', 'performer': [PATH], 'valueQuantity': qty(value, unit),
         'referenceRange': [rr]}
    if interp:
        r['interpretation'] = [interp]
    add(r, [USCORE + 'us-core-observation-lab'])

write_bundle(OUT, 'mtb-l2345-transaction')
