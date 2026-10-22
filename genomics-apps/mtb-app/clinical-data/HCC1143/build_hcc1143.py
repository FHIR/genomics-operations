"""Generate the mCODE transaction bundle for Betty Breast (MRN HCC1143).

Usage: python build_hcc1143.py [output_path]
Writes HCC1143-transaction-bundle.json next to this script by default. The clinical
story behind the data is in HCC1143-scenario.md; edit this script, not the JSON.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # shared mcode_helpers.py

from mcode_helpers import (  # noqa: E402
    ICD10, ICDO3, LOINC, MCODE, RIGHT, RXN, SCT, USCORE, add, cc, condition_base, ext,
    laterality, obs_cat, qty, specimen_type, write_bundle,
)

OUT = sys.argv[1] if len(sys.argv) > 1 else Path(__file__).with_name('HCC1143-transaction-bundle.json')

PID = 'mtb-hcc1143'
PAT = {'reference': f'Patient/{PID}', 'display': 'Betty Breast'}
ONC = {'reference': 'Practitioner/mtb-hcc1143-oncologist', 'display': 'Dr. Olivia Oncologist'}
SURG = {'reference': 'Practitioner/mtb-hcc1143-surgeon', 'display': 'Dr. Sam Surgeon'}
PATH = {'reference': 'Practitioner/mtb-hcc1143-pathologist', 'display': 'Dr. Paul Pathologist'}
RAD = {'reference': 'Practitioner/mtb-hcc1143-radiologist', 'display': 'Dr. Rita Radiologist'}
RADONC = {'reference': 'Practitioner/mtb-hcc1143-radiation-oncologist', 'display': 'Dr. Ray Radonc'}
PRIMARY = {'reference': 'Condition/mtb-hcc1143-primary-cancer', 'display': 'Infiltrating duct carcinoma of right breast'}

# ---------------------------------------------------------------- Patient & practitioners

add({
    'resourceType': 'Patient', 'id': PID,
    'extension': [{'url': 'http://hl7.org/fhir/us/core/StructureDefinition/us-core-birthsex', 'valueCode': 'F'}],
    'identifier': [{
        'use': 'usual',
        'type': cc('http://terminology.hl7.org/CodeSystem/v2-0203', 'MR', 'Medical Record Number'),
        'system': 'http://elimu.io/fhir/mrn',
        'value': 'HCC1143',
    }],
    'active': True,
    'name': [{'use': 'official', 'family': 'Breast', 'given': ['Betty']}],
    'gender': 'female',
    'birthDate': '1972-01-18',
}, [MCODE + 'mcode-cancer-patient'])

for pid, family, given, specialty_code, specialty in [
    ('oncologist', 'Oncologist', 'Olivia', '394593009', 'Medical oncology'),
    ('surgeon', 'Surgeon', 'Sam', '394609007', 'General surgery'),
    ('pathologist', 'Pathologist', 'Paul', '394915009', 'General pathology (specialty)'),
    ('radiologist', 'Radiologist', 'Rita', '394914008', 'Radiology - specialty'),
    ('radiation-oncologist', 'Radonc', 'Ray', '419815003', 'Radiation oncology'),
]:
    add({
        'resourceType': 'Practitioner', 'id': f'mtb-hcc1143-{pid}',
        'identifier': [{'system': 'http://elimu.io/fhir/practitioner-id', 'value': f'mtb-hcc1143-{pid}'}],
        'name': [{'family': family, 'given': [given], 'prefix': ['Dr.']}],
        'qualification': [{'code': cc(SCT, specialty_code, specialty)}],
    }, [USCORE + 'us-core-practitioner'])

# ---------------------------------------------------------------- Specimens

SPEC_CORE = add({
    'resourceType': 'Specimen', 'id': 'mtb-hcc1143-specimen-breast-core',
    'status': 'unavailable',
    'type': specimen_type('TISS', 'Tissue', '258415003', 'Biopsy specimen'),
    'note': [{'text': 'Right breast 10 o\'clock mass, ultrasound-guided core needle biopsy'}],
    'subject': PAT,
    'collection': {
        'collector': RAD, 'collectedDateTime': '2024-04-26',
        'method': cc(SCT, '432550005', 'Core needle biopsy of breast using ultrasound guidance'),
        'bodySite': cc(SCT, '110496004', 'Upper outer quadrant of right breast'),
    },
}, [MCODE + 'mcode-human-specimen'])

SPEC_LUMP = add({
    'resourceType': 'Specimen', 'id': 'mtb-hcc1143-specimen-partial-mastectomy',
    'status': 'available',
    'type': specimen_type('TISS', 'Tissue', '119376003', 'Tissue specimen'),
    'note': [{'text': 'Right partial mastectomy specimen; source of the tumor genomic profiling'}],
    'subject': PAT,
    'collection': {
        'collector': SURG, 'collectedDateTime': '2024-05-20',
        'method': cc(SCT, '736754007', 'Excision of part of right breast'),
        'bodySite': cc(SCT, '110496004', 'Upper outer quadrant of right breast'),
    },
}, [MCODE + 'mcode-human-specimen'])

SPEC_LUNG = add({
    'resourceType': 'Specimen', 'id': 'mtb-hcc1143-specimen-lung-core',
    'status': 'available',
    'type': specimen_type('TISS', 'Tissue', '258415003', 'Biopsy specimen'),
    'note': [{'text': 'Right upper lobe nodule, CT-guided core biopsy'}],
    'subject': PAT,
    'collection': {
        'collector': RAD, 'collectedDateTime': '2026-09-09',
        'method': cc(SCT, '6491000087100', 'CT guided biopsy of right lung'),
        'bodySite': cc(SCT, '3341006', 'Right lung structure'),
    },
}, [MCODE + 'mcode-human-specimen'])

# ---------------------------------------------------------------- Diagnosis

AJCC8 = cc(SCT, '897275008', 'American Joint Committee on Cancer, Cancer Staging Manual, 8th edition neoplasm staging system')

add({
    'resourceType': 'Condition', 'id': 'mtb-hcc1143-primary-cancer',
    'extension': [ext('mcode-histology-morphology-behavior', 'valueCodeableConcept', {'coding': [
        {'system': SCT, 'code': '82711006', 'display': 'Infiltrating duct carcinoma'},
        {'system': ICDO3, 'code': '8500/3', 'display': 'Infiltrating duct carcinoma, NOS'},
    ]})],
    **condition_base('recurrence'),
    'code': {'coding': [
        {'system': SCT, 'code': '408643008', 'display': 'Infiltrating duct carcinoma of breast'},  # in the VSAC Breast Cancer value set
        {'system': ICD10, 'code': 'C50.411', 'display': 'Malignant neoplasm of upper-outer quadrant of right female breast'},
    ], 'text': 'Infiltrating duct carcinoma of right breast'},
    'bodySite': [{
        'extension': laterality(*RIGHT),
        'coding': [
            {'system': SCT, 'code': '110496004', 'display': 'Upper outer quadrant of right breast'},
            {'system': ICDO3, 'code': 'C50.4', 'display': 'Upper-outer quadrant of breast'},
        ],
    }],
    'subject': PAT,
    'onsetDateTime': '2024-04-26',
    'recordedDate': '2024-04-26',
    'asserter': ONC,
    'stage': [
        {
            'summary': cc(SCT, '1222766008', 'IIA (AJCC)', 'Stage IIA (AJCC)'),
            'assessment': [{'reference': 'Observation/mtb-hcc1143-tnm-pathologic-stage-group'}],
            'type': AJCC8,
        },
        {
            'summary': cc(SCT, '1222724007', 'IA (AJCC)', 'Stage IA (AJCC)'),
            'assessment': [{'reference': 'Observation/mtb-hcc1143-tnm-clinical-stage-group'}],
            'type': AJCC8,
        },
    ],
    'evidence': [{'detail': [
        {'reference': 'DiagnosticReport/mtb-hcc1143-path-breast-core'},
        {'reference': 'DiagnosticReport/mtb-hcc1143-path-partial-mastectomy'},
    ]}],
    'note': [{'text': 'Triple-negative (ER 0%, PR 0%, HER2 IHC 0). Recurrence with lung and mediastinal nodal metastases, 2026-08.'}],
}, [MCODE + 'mcode-primary-cancer-condition'])

related_primary = {'url': 'http://hl7.org/fhir/StructureDefinition/condition-related', 'valueReference': PRIMARY}

add({
    'resourceType': 'Condition', 'id': 'mtb-hcc1143-lung-metastases',
    'extension': [related_primary],
    **condition_base(),
    'code': {'coding': [
        {'system': SCT, 'code': '94391008', 'display': 'Metastatic malignant neoplasm to lung'},
        {'system': ICD10, 'code': 'C78.00', 'display': 'Secondary malignant neoplasm of unspecified lung'},
    ]},
    'bodySite': [cc(SCT, '3341006', 'Right lung structure'), cc(SCT, '44029006', 'Left lung structure')],
    'subject': PAT,
    'onsetDateTime': '2026-08-24',
    'recordedDate': '2026-09-12',
    'asserter': ONC,
    'evidence': [{'detail': [
        {'reference': 'DiagnosticReport/mtb-hcc1143-ct-chest-2026-08-24'},
        {'reference': 'DiagnosticReport/mtb-hcc1143-path-lung-biopsy'},
    ]}],
    'note': [{'text': 'Bilateral pulmonary nodules, largest 2.1 cm right upper lobe (biopsy-proven).'}],
}, [MCODE + 'mcode-secondary-cancer-condition'])

add({
    'resourceType': 'Condition', 'id': 'mtb-hcc1143-mediastinal-node-metastases',
    'extension': [related_primary],
    **condition_base(),
    'code': {'coding': [
        {'system': SCT, 'code': '94351005', 'display': 'Metastatic malignant neoplasm to intrathoracic lymph nodes'},
        {'system': ICD10, 'code': 'C77.1', 'display': 'Secondary and unspecified malignant neoplasm of intrathoracic lymph nodes'},
    ]},
    'subject': PAT,
    'onsetDateTime': '2026-08-24',
    'recordedDate': '2026-09-12',
    'asserter': ONC,
    'evidence': [{'detail': [{'reference': 'DiagnosticReport/mtb-hcc1143-pet-ct-2026-08-31'}]}],
    'note': [{'text': 'Subcarinal and right paratracheal nodes, FDG-avid.'}],
}, [MCODE + 'mcode-secondary-cancer-condition'])

# ---------------------------------------------------------------- Staging (AJCC 8th, anatomic): clinical IA, pathologic IIA

def tnm(rid, date, code, value, profile, note=None):
    r = {'resourceType': 'Observation', 'id': rid, 'status': 'final', 'category': obs_cat('exam', 'Exam'),
         'code': code, 'subject': PAT, 'focus': [PRIMARY], 'effectiveDateTime': date, 'performer': [ONC],
         'valueCodeableConcept': value, 'method': AJCC8}
    if note:
        r['note'] = [{'text': note}]
    return add(r, [MCODE + profile])


c_members = [
    tnm('mtb-hcc1143-tnm-clinical-t', '2024-04-30', cc(SCT, '399504009', 'cT category'),
        cc(SCT, '1228899006', 'cT1c (AJCC)'), 'mcode-tnm-primary-tumor-category', '1.9 cm on ultrasound'),
    tnm('mtb-hcc1143-tnm-clinical-n', '2024-04-30', cc(SCT, '399534004', 'cN category'),
        cc(SCT, '1229967007', 'cN0 (AJCC)'), 'mcode-tnm-regional-nodes-category', 'No suspicious axillary nodes on ultrasound'),
    tnm('mtb-hcc1143-tnm-clinical-m', '2024-04-30', cc(SCT, '399387003', 'cM category'),
        cc(SCT, '1229901006', 'cM0 (AJCC)'), 'mcode-tnm-distant-metastases-category'),
]
add({'resourceType': 'Observation', 'id': 'mtb-hcc1143-tnm-clinical-stage-group', 'status': 'final',
     'category': obs_cat('exam', 'Exam'), 'code': cc(SCT, '399537006', 'Clinical TNM stage grouping'),
     'subject': PAT, 'focus': [PRIMARY], 'effectiveDateTime': '2024-04-30', 'performer': [ONC],
     'valueCodeableConcept': cc(SCT, '1222724007', 'IA (AJCC)'), 'method': AJCC8, 'hasMember': c_members},
    [MCODE + 'mcode-tnm-stage-group'])

p_members = [
    tnm('mtb-hcc1143-tnm-pathologic-t', '2024-05-28', cc(SCT, '384625004', 'pT category'),
        cc(SCT, '1229852009', 'pT2 (AJCC)'), 'mcode-tnm-primary-tumor-category', '2.3 cm invasive carcinoma'),
    tnm('mtb-hcc1143-tnm-pathologic-n', '2024-05-28', cc(SCT, '371494008', 'pN category'),
        cc(SCT, '1229947003', 'pN0 (AJCC)'), 'mcode-tnm-regional-nodes-category', '0 of 3 sentinel lymph nodes involved'),
    c_members[2],
]
add({'resourceType': 'Observation', 'id': 'mtb-hcc1143-tnm-pathologic-stage-group', 'status': 'final',
     'category': obs_cat('exam', 'Exam'), 'code': cc(SCT, '399588009', 'Pathologic TNM stage grouping'),
     'subject': PAT, 'focus': [PRIMARY], 'effectiveDateTime': '2024-05-28', 'performer': [ONC],
     'valueCodeableConcept': cc(SCT, '1222766008', 'IIA (AJCC)'), 'method': AJCC8, 'hasMember': p_members,
     'note': [{'text': 'Anatomic stage group. (AJCC 8th edition prognostic stage for a grade 3 triple-negative T2N0 tumor is higher.)'}]},
    [MCODE + 'mcode-tnm-stage-group'])

add({'resourceType': 'Observation', 'id': 'mtb-hcc1143-histologic-grade', 'status': 'final',
     'category': obs_cat('laboratory', 'Laboratory'),
     'code': {'coding': [{'system': 'http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl', 'code': 'C18000',
                          'display': 'Histologic Grade'}]},
     'subject': PAT, 'focus': [PRIMARY], 'effectiveDateTime': '2024-05-20', 'performer': [PATH],
     'valueCodeableConcept': cc(SCT, '1155707008', 'High histologic grade'),
     'method': cc(SCT, '449205006', 'Nottingham histologic grading system'),
     'specimen': SPEC_LUMP, 'note': [{'text': 'Nottingham grade 3 (tubules 3, pleomorphism 3, mitoses 3; score 9/9)'}]},
    [MCODE + 'mcode-histologic-grade'])

# ---------------------------------------------------------------- Tumor and tumor size

TUMOR = add({
    'resourceType': 'BodyStructure', 'id': 'mtb-hcc1143-tumor-breast',
    'identifier': [{'use': 'usual', 'type': cc('http://hl7.org/fhir/resource-types', 'BodyStructure', 'BodyStructure'),
                    'system': 'http://elimu.io/fhir/tumor-id', 'value': 'HCC1143-right-breast-primary'}],
    'morphology': cc(SCT, '82711006', 'Infiltrating duct carcinoma'),
    'location': cc(SCT, '110496004', 'Upper outer quadrant of right breast'),
    'locationQualifier': [cc(SCT, '24028007', 'Right')],
    'description': 'Primary tumor, right breast 10 o\'clock',
    'patient': PAT,
}, [MCODE + 'mcode-tumor'])


def tumor_size(rid, date, longest, other, method, performer, specimen=None):
    r = {'resourceType': 'Observation', 'id': rid, 'status': 'final',
         'category': obs_cat('laboratory' if specimen else 'imaging', 'Laboratory' if specimen else 'Imaging'),
         'code': cc(LOINC, '21889-1', 'Size Tumor'),
         'subject': PAT, 'effectiveDateTime': date, 'performer': [performer], 'method': method,
         'component': [
             {'code': cc(LOINC, '33728-7', 'Size.maximum dimension in Tumor'), 'valueQuantity': qty(longest, 'cm')},
             {'code': cc(LOINC, '33729-5', 'Size additional dimension in Tumor'), 'valueQuantity': qty(other, 'cm')},
         ]}
    if specimen:
        r['specimen'] = specimen
    else:
        r['focus'] = [TUMOR]
    add(r, [MCODE + 'mcode-tumor-size'])


tumor_size('mtb-hcc1143-tumor-size-2024-04-22', '2024-04-22', 1.9, 1.5, cc(SCT, '16310003', 'Ultrasonography'), RAD)
tumor_size('mtb-hcc1143-tumor-size-2024-05-20', '2024-05-20', 2.3, 1.8,
           cc(SCT, '787377000', 'Gross examination and sampling of tissue specimen'), PATH, SPEC_LUMP)

# ---------------------------------------------------------------- Reports

DR_CAT_RAD = cc(LOINC, 'LP29684-5', 'Radiology')
DR_CAT_PATH = cc(LOINC, 'LP7839-6', 'Pathology')


def report(rid, date, cat, code, conclusion, performer, specimen=None, conclusion_code=None):
    r = {'resourceType': 'DiagnosticReport', 'id': rid, 'status': 'final', 'category': [cat], 'code': code,
         'subject': PAT, 'effectiveDateTime': date, 'issued': f'{date}T17:00:00-04:00',
         'performer': [performer], 'conclusion': conclusion}
    if specimen:
        r['specimen'] = [specimen]
    if conclusion_code:
        r['conclusionCode'] = [conclusion_code]
    add(r, [USCORE + 'us-core-diagnosticreport-note'])


IDC = cc(SCT, '82711006', 'Infiltrating duct carcinoma')

report('mtb-hcc1143-mammogram-2024-04-22', '2024-04-22', DR_CAT_RAD, cc(LOINC, '24605-8', 'MG Breast Diagnostic'),
       'Irregular, spiculated mass at 10 o\'clock in the right breast, posterior depth, with associated architectural '
       'distortion. No suspicious calcifications. Left breast unremarkable. BI-RADS 5.', RAD)
report('mtb-hcc1143-us-breast-2024-04-22', '2024-04-22', DR_CAT_RAD, cc(LOINC, '24601-7', 'US Breast'),
       'Right breast 10 o\'clock, 6 cm from the nipple: irregular hypoechoic mass with angular margins, 1.9 x 1.5 cm. '
       'Right axilla: no morphologically abnormal lymph nodes. BI-RADS 5; ultrasound-guided biopsy recommended.', RAD)
report('mtb-hcc1143-path-breast-core', '2024-04-26', DR_CAT_PATH, cc(LOINC, '11526-1', 'Pathology study'),
       'Right breast 10 o\'clock, ultrasound-guided core biopsy: invasive ductal carcinoma, Nottingham grade 3. '
       'ER negative (0%), PR negative (0%), HER2 negative (IHC 0). Ki-67 approximately 70%.', PATH, SPEC_CORE, IDC)
report('mtb-hcc1143-path-partial-mastectomy', '2024-05-20', DR_CAT_PATH, cc(LOINC, '11526-1', 'Pathology study'),
       'Right breast partial mastectomy and sentinel lymph node biopsy: invasive ductal carcinoma, Nottingham grade 3 '
       '(score 9/9), 2.3 cm. No lymphovascular invasion. Margins negative (closest 4 mm, posterior). '
       'Sentinel lymph nodes: 0 of 3 involved. pT2 pN0.', PATH, SPEC_LUMP, IDC)
report('mtb-hcc1143-germline-2024-06-07', '2024-06-07', DR_CAT_PATH, cc(LOINC, '51969-4', 'Genetic analysis report'),
       'Germline multigene hereditary cancer panel (peripheral blood): no pathogenic or likely pathogenic variants '
       'identified, including BRCA1, BRCA2, PALB2, TP53, PTEN, CHEK2 and ATM.', PATH)
report('mtb-hcc1143-molecular-2024-06-14', '2024-06-14', DR_CAT_PATH, cc(LOINC, '51969-4', 'Genetic analysis report'),
       'Tumor genomic profiling (partial mastectomy specimen, with matched normal): TP53 c.743G>A (p.Arg248Gln), '
       'homozygous. No pathogenic alterations in PIK3CA, AKT1, PTEN, BRCA1 or BRCA2. Microsatellite stable. '
       'Detailed variant-level results are maintained on the genomics server, not in this report.', PATH, SPEC_LUMP)
report('mtb-hcc1143-mammogram-2025-05-19', '2025-05-19', DR_CAT_RAD, cc(LOINC, '24605-8', 'MG Breast Diagnostic'),
       'Post-lumpectomy and radiation changes in the right breast. No suspicious mass or calcifications. '
       'Left breast unremarkable. BI-RADS 2.', RAD)
report('mtb-hcc1143-ct-chest-2026-08-24', '2026-08-24', DR_CAT_RAD, cc(LOINC, '24628-0', 'CT Chest W contrast IV'),
       'Multiple new bilateral pulmonary nodules, the largest 2.1 cm in the right upper lobe, highly suspicious for '
       'metastases. Enlarged subcarinal (1.6 cm) and right paratracheal (1.3 cm) lymph nodes. No pleural effusion. '
       'Post-surgical changes in the right breast.', RAD)
report('mtb-hcc1143-pet-ct-2026-08-31', '2026-08-31', DR_CAT_RAD,
       cc(LOINC, '81554-8', 'PET+CT Guidance for localization of tumor of Skull base to mid-thigh-- W 18F-FDG IV'),
       'FDG-avid bilateral pulmonary nodules (SUVmax 8.9) and subcarinal and right paratracheal lymph nodes '
       '(SUVmax 7.2). No FDG-avid hepatic, osseous or adrenal disease. No residual disease in the right breast.', RAD)
report('mtb-hcc1143-mri-brain-2026-09-02', '2026-09-02', DR_CAT_RAD, cc(LOINC, '24587-8', 'MR Brain WO and W contrast IV'),
       'No evidence of intracranial metastatic disease.', RAD)
report('mtb-hcc1143-path-lung-biopsy', '2026-09-09', DR_CAT_PATH, cc(LOINC, '11526-1', 'Pathology study'),
       'Right upper lobe nodule, CT-guided core biopsy: metastatic carcinoma consistent with the known breast primary '
       '(GATA3 positive, TTF-1 negative). ER negative (0%), PR negative (0%), HER2 negative (IHC 0). '
       'PD-L1 (22C3) combined positive score 12.', PATH, SPEC_LUNG,
       cc(SCT, '82711006', 'Infiltrating duct carcinoma'))

# ---------------------------------------------------------------- Tumor markers

POS_NEG = 'http://loinc.org'
NEG = cc(POS_NEG, 'LA6577-6', 'Negative')


def tumor_marker(rid, date, code, value, specimen, note=None, interpretation=None, ref=None):
    r = {'resourceType': 'Observation', 'id': rid, 'status': 'final',
         'extension': [ext('mcode-related-condition', 'valueReference', PRIMARY)],
         'category': obs_cat('laboratory', 'Laboratory'), 'code': code, 'subject': PAT,
         'effectiveDateTime': date, 'performer': [PATH], **value}
    if specimen:
        r['specimen'] = specimen
    if note:
        r['note'] = [{'text': note}]
    if interpretation:
        r['interpretation'] = [interpretation]
    if ref:
        r['referenceRange'] = [ref]
    add(r, [MCODE + 'mcode-tumor-marker-test'])


ER = cc(LOINC, '85337-4', 'Estrogen receptor Ag [Presence] in Breast cancer specimen by Immune stain')
PR = cc(LOINC, '85339-0', 'Progesterone receptor Ag [Presence] in Breast cancer specimen by Immune stain')
HER2 = cc(LOINC, '85319-2', 'HER2 Ag [Presence] in Breast cancer specimen by Immune stain')

for suffix, date, spec in [('2024-04-26', '2024-04-26', SPEC_CORE), ('2026-09-09', '2026-09-09', SPEC_LUNG)]:
    tumor_marker(f'mtb-hcc1143-er-{suffix}', date, ER, {'valueCodeableConcept': NEG}, spec, 'ER 0% of tumor nuclei staining')
    tumor_marker(f'mtb-hcc1143-pr-{suffix}', date, PR, {'valueCodeableConcept': NEG}, spec, 'PR 0% of tumor nuclei staining')
    tumor_marker(f'mtb-hcc1143-her2-{suffix}', date, HER2, {'valueCodeableConcept': NEG}, spec, 'HER2 IHC score 0')

tumor_marker('mtb-hcc1143-ki67-2024-04-26', '2024-04-26',
             cc(LOINC, '85330-9', 'Cells.Ki-67 nuclear Ag/cells in Breast cancer specimen by Immune stain'),
             {'valueQuantity': qty(70, '%')}, SPEC_CORE)
tumor_marker('mtb-hcc1143-pdl1-cps-2026-09-09', '2026-09-09',
             cc(LOINC, '85147-7', 'PD-L1 by clone 22C3 in Tissue by Immune stain Report'),
             {'valueString': 'Combined positive score (CPS) 12'}, SPEC_LUNG,
             'PD-L1 IHC 22C3 pharmDx on lung metastasis. CPS >= 10.')

HIGH = cc('http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation', 'H', 'High')
for date, value in [('2024-04-26', 24), ('2025-05-19', 19), ('2026-08-28', 71)]:
    tumor_marker(f'mtb-hcc1143-ca15-3-{date}', date,
                 cc(LOINC, '6875-9', 'Cancer Ag 15-3 [Units/volume] in Serum or Plasma'),
                 {'valueQuantity': qty(value, 'U/mL', '[iU]/mL')}, None,
                 interpretation=HIGH if value > 30 else None,
                 ref={'high': qty(30, 'U/mL', '[iU]/mL'), 'text': '<=30 U/mL'})

# ---------------------------------------------------------------- Disease status

add({'resourceType': 'Observation', 'id': 'mtb-hcc1143-disease-status-2026-08-24',
     'extension': [ext('mcode-cancer-disease-status-evidence-type', 'valueCodeableConcept',
                       cc(SCT, '363679005', 'Imaging'))],
     'status': 'final', 'category': obs_cat('exam', 'Exam'),
     'code': cc(LOINC, '97509-4', 'Cancer disease progression'),
     'subject': PAT, 'focus': [PRIMARY], 'effectiveDateTime': '2026-08-24', 'performer': [ONC],
     'valueCodeableConcept': cc(SCT, '271299001', "Patient's condition worsened"),
     'note': [{'text': 'Distant recurrence: bilateral lung and mediastinal nodal metastases.'}]},
    [MCODE + 'mcode-cancer-disease-status'])

# ---------------------------------------------------------------- ECOG

for date, score, la, display in [
    ('2024-04-26', 0, 'LA9622-7', 'Fully active, able to carry on all pre-disease performance without restriction'),
    ('2026-09-15', 1, 'LA9623-5', 'Restricted in physically strenuous activity but ambulatory and able to carry out work of a light or sedentary nature, e.g., light house work, office work'),
]:
    add({'resourceType': 'Observation', 'id': f'mtb-hcc1143-ecog-{date}', 'status': 'final',
         'category': obs_cat('survey', 'Survey'),
         'code': cc(LOINC, '89247-1', 'ECOG Performance Status score'),
         'subject': PAT, 'effectiveDateTime': date, 'performer': [ONC], 'valueInteger': score,
         'interpretation': [cc(LOINC, la, display)]},
        [MCODE + 'mcode-ecog-performance-status'])

# ---------------------------------------------------------------- Cancer-related medications (adjuvant)

CURATIVE = ext('mcode-procedure-intent', 'valueCodeableConcept', cc(SCT, '373808002', 'Curative - procedure intent'))
SUPPORTIVE = ext('mcode-procedure-intent', 'valueCodeableConcept', cc(SCT, '399707004', 'Supportive - procedure intent'))
IV = cc(SCT, '47625008', 'Intravenous route')
SUBQ = cc(SCT, '34206005', 'Subcutaneous route')


def med_request(rid, rxcui, display, start, end, dosage_text, timing, route, intent_ext=CURATIVE, dose=None, note=None):
    dosage = {'text': dosage_text, 'route': route,
              'timing': {'repeat': {'boundsPeriod': {'start': start, 'end': end}, **timing}}}
    if dose:
        dosage['doseAndRate'] = [dose]
    r = {'resourceType': 'MedicationRequest', 'id': rid, 'extension': [intent_ext], 'status': 'completed',
         'intent': 'order', 'category': [cc('http://terminology.hl7.org/CodeSystem/medicationrequest-category',
                                            'outpatient', 'Outpatient')],
         'medicationCodeableConcept': cc(RXN, rxcui, display), 'subject': PAT, 'authoredOn': start,
         'requester': ONC, 'reasonReference': [PRIMARY], 'dosageInstruction': [dosage]}
    if note:
        r['note'] = [{'text': note}]
    add(r, [MCODE + 'mcode-cancer-related-medication-request'])


DD = {'count': 4, 'frequency': 1, 'period': 14, 'periodUnit': 'd'}
med_request('mtb-hcc1143-medrx-doxorubicin', '1790099', '10 ML doxorubicin hydrochloride 2 MG/ML Injection',
            '2024-06-24', '2024-08-05', 'Doxorubicin 60 mg/m2 IV on day 1 every 14 days x 4 cycles (dose-dense AC)', DD, IV,
            dose={'doseQuantity': qty(60, 'mg/m2', 'mg/m2')},
            note='Adjuvant dose-dense AC, cycles 2024-06-24, 07-08, 07-22, 08-05.')
med_request('mtb-hcc1143-medrx-cyclophosphamide', '1734917', 'cyclophosphamide 500 MG Injection',
            '2024-06-24', '2024-08-05', 'Cyclophosphamide 600 mg/m2 IV on day 1 every 14 days x 4 cycles (dose-dense AC)', DD, IV,
            dose={'doseQuantity': qty(600, 'mg/m2', 'mg/m2')})
med_request('mtb-hcc1143-medrx-pegfilgrastim', '727539', '0.6 ML pegfilgrastim 10 MG/ML Prefilled Syringe',
            '2024-06-25', '2024-08-06', 'Pegfilgrastim 6 mg subcutaneously on day 2 of each dose-dense AC cycle',
            DD, SUBQ, intent_ext=SUPPORTIVE)
med_request('mtb-hcc1143-medrx-paclitaxel', '312199', 'paclitaxel 6 MG/ML Injectable Solution',
            '2024-08-19', '2024-11-04', 'Paclitaxel 80 mg/m2 IV weekly x 12',
            {'count': 12, 'frequency': 1, 'period': 1, 'periodUnit': 'wk'}, IV,
            dose={'doseQuantity': qty(80, 'mg/m2', 'mg/m2')}, note='Adjuvant weekly paclitaxel after dose-dense AC.')

# ---------------------------------------------------------------- Radiotherapy (adjuvant whole breast)

RT_VOLUME = add({
    'resourceType': 'BodyStructure', 'id': 'mtb-hcc1143-rt-volume-right-breast',
    'identifier': [{'use': 'usual', 'system': 'http://elimu.io/fhir/rt-volume', 'value': 'PTV_R_Breast'}],
    'morphology': cc(SCT, '228793007', 'Planning target volume'),
    'location': cc(SCT, '76752008', 'Breast structure'),
    'locationQualifier': [cc(SCT, '24028007', 'Right')],
    'description': 'Right whole breast',
    'patient': PAT,
}, [MCODE + 'mcode-radiotherapy-volume'])

add({
    'resourceType': 'Procedure', 'id': 'mtb-hcc1143-rt-course-right-breast',
    'extension': [
        CURATIVE,
        {'url': MCODE + 'mcode-radiotherapy-modality-and-technique', 'extension': [
            {'url': MCODE + 'mcode-radiotherapy-modality', 'valueCodeableConcept':
                cc(SCT, '1156506007', 'External beam radiation therapy using photons')},
            {'url': MCODE + 'mcode-radiotherapy-technique', 'valueCodeableConcept':
                cc(SCT, '1162782007', 'Three dimensional external beam radiation therapy')},
        ]},
        {'url': MCODE + 'mcode-radiotherapy-sessions', 'valueUnsignedInt': 15},
        {'url': MCODE + 'mcode-radiotherapy-dose-delivered-to-volume', 'extension': [
            {'url': 'volume', 'valueReference': RT_VOLUME},
            {'url': 'totalDoseDelivered', 'valueQuantity': qty(4005, 'cGy')},
            {'url': 'fractionsDelivered', 'valueUnsignedInt': 15},
        ]},
    ],
    'status': 'completed',
    'code': cc(SCT, '1217123003', 'Radiotherapy course of treatment'),
    'subject': PAT,
    'performedPeriod': {'start': '2024-11-25', 'end': '2024-12-13'},
    'performer': [{'actor': RADONC}],
    'reasonReference': [PRIMARY],
    'bodySite': [{'extension': laterality(*RIGHT), **cc(SCT, '76752008', 'Breast structure')}],
    'note': [{'text': 'Adjuvant hypofractionated whole-breast radiotherapy, 40.05 Gy in 15 fractions.'}],
}, [MCODE + 'mcode-radiotherapy-course-summary'])

# ---------------------------------------------------------------- Procedures (biopsies, surgery)


def procedure(rid, date, code, site, report_id, performer, profile, intent=None):
    r = {'resourceType': 'Procedure', 'id': rid, 'status': 'completed', 'code': code, 'subject': PAT,
         'performedDateTime': date, 'performer': [{'actor': performer}], 'reasonReference': [PRIMARY],
         'bodySite': [site], 'report': [{'reference': f'DiagnosticReport/{report_id}'}]}
    if intent:
        r['extension'] = [intent]
    add(r, [profile])


procedure('mtb-hcc1143-proc-breast-core', '2024-04-26',
          cc(SCT, '432550005', 'Core needle biopsy of breast using ultrasound guidance'),
          cc(SCT, '110496004', 'Upper outer quadrant of right breast'), 'mtb-hcc1143-path-breast-core', RAD,
          USCORE + 'us-core-procedure')
procedure('mtb-hcc1143-proc-partial-mastectomy', '2024-05-20', cc(SCT, '736754007', 'Excision of part of right breast'),
          cc(SCT, '110496004', 'Upper outer quadrant of right breast'), 'mtb-hcc1143-path-partial-mastectomy', SURG,
          MCODE + 'mcode-cancer-related-surgical-procedure', CURATIVE)
procedure('mtb-hcc1143-proc-sentinel-node', '2024-05-20', cc(SCT, '396487001', 'Sentinel lymph node biopsy'),
          {'extension': laterality(*RIGHT), **cc(SCT, '68171009', 'Axillary lymph node structure')},
          'mtb-hcc1143-path-partial-mastectomy', SURG, MCODE + 'mcode-cancer-related-surgical-procedure', CURATIVE)
procedure('mtb-hcc1143-proc-lung-biopsy', '2026-09-09', cc(SCT, '6491000087100', 'CT guided biopsy of right lung'),
          cc(SCT, '3341006', 'Right lung structure'), 'mtb-hcc1143-path-lung-biopsy', RAD, USCORE + 'us-core-procedure')

# ---------------------------------------------------------------- Other conditions, medications, allergies, history

hypothyroid = add({'resourceType': 'Condition', 'id': 'mtb-hcc1143-hypothyroidism', **condition_base(),
                   'code': cc(SCT, '40930008', 'Hypothyroidism'),
                   'subject': PAT, 'onsetDateTime': '2015', 'recordedDate': '2015-03-02'},
                  [USCORE + 'us-core-condition-problems-health-concerns'])

add({'resourceType': 'MedicationRequest', 'id': 'mtb-hcc1143-medrx-levothyroxine', 'status': 'active', 'intent': 'order',
     'medicationCodeableConcept': cc(RXN, '966222', 'levothyroxine sodium 0.075 MG Oral Tablet'),
     'subject': PAT, 'authoredOn': '2026-01-12', 'requester': ONC, 'reasonReference': [hypothyroid],
     'dosageInstruction': [{'text': 'Levothyroxine 75 mcg by mouth once daily', 'route': cc(SCT, '26643006', 'Oral route'),
                            'timing': {'repeat': {'frequency': 1, 'period': 1, 'periodUnit': 'd'}}}]},
    [USCORE + 'us-core-medicationrequest'])

add({'resourceType': 'AllergyIntolerance', 'id': 'mtb-hcc1143-no-known-allergy',
     'clinicalStatus': cc('http://terminology.hl7.org/CodeSystem/allergyintolerance-clinical', 'active', 'Active'),
     'verificationStatus': cc('http://terminology.hl7.org/CodeSystem/allergyintolerance-verification', 'confirmed', 'Confirmed'),
     'code': cc(SCT, '716186003', 'No known allergy'),
     'patient': PAT, 'recordedDate': '2024-04-26'},
    [USCORE + 'us-core-allergyintolerance'])

add({'resourceType': 'Observation', 'id': 'mtb-hcc1143-smoking-status', 'status': 'final',
     'category': obs_cat('social-history', 'Social History'),
     'code': cc(LOINC, '72166-2', 'Tobacco smoking status'), 'subject': PAT,
     'effectiveDateTime': '2024-04-26', 'performer': [ONC],
     'valueCodeableConcept': cc(SCT, '8517006', 'Ex-smoker'),
     'note': [{'text': 'Quit 2010; about 10 pack-years.'}]},
    [USCORE + 'us-core-smokingstatus'])

add({'resourceType': 'FamilyMemberHistory', 'id': 'mtb-hcc1143-family-history-mother', 'status': 'completed',
     'patient': PAT, 'date': '2024-04-26',
     'relationship': cc('http://terminology.hl7.org/CodeSystem/v3-RoleCode', 'MTH', 'mother'),
     'sex': cc('http://hl7.org/fhir/administrative-gender', 'female', 'Female'),
     'condition': [{'code': cc(SCT, '254837009', 'Malignant neoplasm of breast'),
                    'onsetAge': {'value': 64, 'unit': 'yr', 'system': 'http://unitsofmeasure.org', 'code': 'a'}}]})

# ---------------------------------------------------------------- Vital signs

VS_CAT = obs_cat('vital-signs', 'Vital Signs')
for date, weight in [('2024-04-26', 84.0), ('2026-09-15', 80.6)]:
    add({'resourceType': 'Observation', 'id': f'mtb-hcc1143-weight-{date}', 'status': 'final', 'category': VS_CAT,
         'code': cc(LOINC, '29463-7', 'Body weight'), 'subject': PAT, 'effectiveDateTime': date, 'performer': [ONC],
         'valueQuantity': qty(weight, 'kg')},
        ['http://hl7.org/fhir/StructureDefinition/bodyweight', USCORE + 'us-core-body-weight'])
add({'resourceType': 'Observation', 'id': 'mtb-hcc1143-height-2024-04-26', 'status': 'final', 'category': VS_CAT,
     'code': cc(LOINC, '8302-2', 'Body height'), 'subject': PAT, 'effectiveDateTime': '2024-04-26', 'performer': [ONC],
     'valueQuantity': qty(165, 'cm')},
    ['http://hl7.org/fhir/StructureDefinition/bodyheight', USCORE + 'us-core-body-height'])

# ---------------------------------------------------------------- Current labs (2026-09-15)

LOW = cc('http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation', 'L', 'Low')
for code, display, value, unit, low, high, interp in [
    ('2160-0', 'Creatinine [Mass/volume] in Serum or Plasma', 0.8, 'mg/dL', 0.5, 1.1, None),
    ('98979-8', 'Glomerular filtration rate [Volume Rate/Area] in Serum, Plasma or Blood by Creatinine-based formula (CKD-EPI 2021)/1.73 sq M', 86, 'mL/min/{1.73_m2}', 60, None, None),
    ('1920-8', 'Aspartate aminotransferase [Enzymatic activity/volume] in Serum or Plasma', 24, 'U/L', 10, 40, None),
    ('1742-6', 'Alanine aminotransferase [Enzymatic activity/volume] in Serum or Plasma', 21, 'U/L', 7, 56, None),
    ('6768-6', 'Alkaline phosphatase [Enzymatic activity/volume] in Serum or Plasma', 88, 'U/L', 44, 147, None),
    ('1975-2', 'Bilirubin.total [Mass/volume] in Serum or Plasma', 0.6, 'mg/dL', 0.1, 1.2, None),
    ('718-7', 'Hemoglobin [Mass/volume] in Blood', 11.6, 'g/dL', 12.0, 15.5, LOW),
    ('6690-2', 'Leukocytes [#/volume] in Blood by Automated count', 6.4, '10*3/uL', 4.0, 11.0, None),
    ('751-8', 'Neutrophils [#/volume] in Blood by Automated count', 4.1, '10*3/uL', 1.8, 7.7, None),
    ('777-3', 'Platelets [#/volume] in Blood by Automated count', 268, '10*3/uL', 150, 400, None),
]:
    rr = {**({'low': qty(low, unit)} if low is not None else {}), **({'high': qty(high, unit)} if high is not None else {})}
    r = {'resourceType': 'Observation', 'id': f'mtb-hcc1143-lab-{code}-2026-09-15', 'status': 'final',
         'category': obs_cat('laboratory', 'Laboratory'), 'code': cc(LOINC, code, display), 'subject': PAT,
         'effectiveDateTime': '2026-09-15', 'performer': [PATH], 'valueQuantity': qty(value, unit),
         'referenceRange': [rr]}
    if interp:
        r['interpretation'] = [interp]
    add(r, [USCORE + 'us-core-observation-lab'])

write_bundle(OUT, 'mtb-hcc1143-transaction')
