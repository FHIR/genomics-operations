"""Build public/data/MTB_TumorTypeCrosswalk.csv from the VSAC value-set spreadsheets in value-sets/.

Each value set maps to a predictor tumorType (blank when the predictor has no matching type) and a
gene-list "Cancer category" (as in public/data/MTB_KB_GeneLists.csv). The spreadsheets are used
exactly as downloaded from VSAC; to refresh, download newer versions and rerun:

    python terminology/build_crosswalk.py
"""
import csv
import warnings
from pathlib import Path

import openpyxl

# VSAC workbooks have no default style; openpyxl warns about it on every load
warnings.filterwarnings('ignore', message='Workbook contains no default style', module='openpyxl')

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'public' / 'data' / 'MTB_TumorTypeCrosswalk.csv'

# value-set file -> (predictor tumorType, gene-list Cancer category)
VALUE_SETS = {
    'BreastCancer.xlsx': ('Breast Cancer', 'Breast Cancer'),
    'NonSmallCellLungCancer.xlsx': ('Non-Small Cell Lung Cancer', 'Non-Small Cell Lung Cancer'),
    'HepatocellularCarcinoma.xlsx': ('', 'Hepatocellular Carcinoma'),
    'ColorectalCancer.xlsx': ('', 'Colorectal Cancer'),
}

# Codes added to a tumor type when its VSAC value set lacks them; their crosswalk rows get valueSetOid "local".
# If a newer VSAC version adds one of these codes, the build fails as a duplicate: remove it from this list.
LOCAL_CODES = [
    # The ASCO colorectal value set has no adenocarcinoma codes, and the TCGA colorectal patients
    # are coded with these three
    ('http://snomed.info/sct', '1701000119104', 'Primary adenocarcinoma of colon (disorder)', 'Colorectal Cancer'),
    ('http://snomed.info/sct', '254582000', 'Adenocarcinoma of rectum (disorder)', 'Colorectal Cancer'),
    ('http://snomed.info/sct', '425178004', 'Adenocarcinoma of rectosigmoid junction (disorder)', 'Colorectal Cancer'),
]

# VSAC code system names -> FHIR system URIs
CODE_SYSTEMS = {
    'SNOMEDCT': 'http://snomed.info/sct',
}


def value_set_oid(workbook):
    for row in workbook['Value Set Info'].iter_rows(values_only=True):
        if row and row[0] == 'OID':
            return row[1]
    raise ValueError('no OID in Value Set Info sheet')


def expansion(workbook):
    rows = workbook['Expansion List'].iter_rows(values_only=True)
    for row in rows:
        if row and row[0] == 'Code':
            header = list(row)
            break
    else:
        raise ValueError('no header row in Expansion List sheet')
    for row in rows:
        if row and row[0]:
            yield dict(zip(header, row))


def main():
    out_rows = []
    seen = {}
    for filename, (predictor_tumor_type, gene_list_category) in VALUE_SETS.items():
        workbook = openpyxl.load_workbook(HERE / 'value-sets' / filename, read_only=True)
        oid = value_set_oid(workbook)
        for entry in expansion(workbook):
            system = CODE_SYSTEMS[str(entry['Code System'])]
            code = str(entry['Code'])
            if (system, code) in seen:
                raise ValueError(f'{code} is in both {seen[(system, code)]} and {filename}')
            seen[(system, code)] = filename
            out_rows.append({
                'system': system,
                'code': code,
                'display': entry['Description'],
                'predictorTumorType': predictor_tumor_type,
                'geneListCategory': gene_list_category,
                'valueSetOid': oid,
            })

    predictor_by_category = {category: predictor for predictor, category in VALUE_SETS.values()}
    for system, code, display, gene_list_category in LOCAL_CODES:
        if (system, code) in seen:
            raise ValueError(f'local code {code} is already in {seen[(system, code)]}')
        seen[(system, code)] = 'LOCAL_CODES'
        out_rows.append({
            'system': system,
            'code': code,
            'display': display,
            'predictorTumorType': predictor_by_category[gene_list_category],
            'geneListCategory': gene_list_category,
            'valueSetOid': 'local',
        })

    with open(OUT, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f'{len(out_rows)} codes -> {OUT}')


if __name__ == '__main__':
    main()
