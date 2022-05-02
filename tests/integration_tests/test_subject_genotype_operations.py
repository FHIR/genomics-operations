import pytest
from tests.utilities import *


"""
Find Subject Variants Tests
---------------------------
"""

def test_find_subject_variants_subject_supplied(client):
    url = find_subject_variants_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_variants_ranges_supplied(client):
    url = find_subject_variants_query('subject=HG00403&testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_variants_1(client):
    url = find_subject_variants_query('subject=HG00403&ranges=NC_000019.10:11089431-11133820&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_find_subject_variants_2(client):
    url = find_subject_variants_query('subject=HG00403&ranges=NC_000019.9:11200138-11244496&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_VARIANTS_OUTPUT_DIR}2.json', response.json)


def test_find_subject_variants_3(client):
    url = find_subject_variants_query('subject=HG00403&includeVariants=true&ranges=NC_000019.10:11089431-11133820')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_VARIANTS_OUTPUT_DIR}3.json', response.json)


def test_find_subject_variants_4(client):
    url = find_subject_variants_query('subject=HG00403&ranges=NC_000019.9:11221327-11221447&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_VARIANTS_OUTPUT_DIR}4.json', response.json)


def test_find_subject_variants_5(client):
    url = find_subject_variants_query('subject=HG00403&ranges=NC_000001.10:69501-69520   ,   NC_000001.10:69501-69521  &testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_VARIANTS_OUTPUT_DIR}5.json', response.json)


"""
Find Subject Specific Variants Tests
------------------------------------
"""

def test_find_subject_specific_variants_subject_supplied(client):
    url = find_subject_variants_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_specific_variants_variants_supplied(client):
    url = find_subject_variants_query('subject=HG00403&testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400

def test_find_subject_specific_variants_failed_normalization(client):
    url = find_subject_variants_query('subject=m123&variants=NM_001385641.1:c.804C>T')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_specific_variants_1(client):
    url = find_subject_specific_variants_query('subject=m123&variants=NC_000005.9:g.112102032A>T')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_SPECIFIC_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_find_subject_specific_variants_2(client):
    url = find_subject_specific_variants_query('subject=m123&variants=NM_001127510.3:c.145A>T')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_SPECIFIC_VARIANTS_OUTPUT_DIR}2.json', response.json)


def test_find_subject_specific_variants_3(client):
    url = find_subject_specific_variants_query('subject=HG00404&variants=NC_000005.9:112154946:CC:C,NC_000005.9:112157612:C:T,NC_000017.10:41246373:G:,NC_000017.10:41243854:AAA:AAAA,NC_000013.10:32900284:C:A,NC_000019.9:11216262:C:A,NC_000014.8:23894101:A:T')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_SPECIFIC_VARIANTS_OUTPUT_DIR}3.json', response.json)


"""
Find Subject Structural Intersecting Variants Tests
---------------------------------------------------
"""

def test_find_subject_structural_intersecting_variants_subject_supplied(client):
    url = find_subject_structural_intersecting_variants_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_structural_intersecting_variants_ranges_supplied(client):
    url = find_subject_structural_intersecting_variants_query('subject=HG00403&testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_structural_intersecting_variants_1(client):
    url = find_subject_structural_intersecting_variants_query('subject=HCC1143&ranges=NC_000007.14:55019016-55211628,NC_000004.12:54727415-54727542&specimenIdentifiers=HCC1143-Sp1')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_find_subject_structural_intersecting_variants_2(client):
    url = find_subject_structural_intersecting_variants_query('subject=HCC1143&ranges=NC_000007.14:55019016-55211628,NC_000004.12:54727415-54727542&specimenIdentifiers=HCC1143-Sp1')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}2.json', response.json)


def test_find_subject_structural_intersecting_variants_3(client):
    url = find_subject_structural_intersecting_variants_query('subject=ABC789&ranges=NC_000002.12:179400709-179483218&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}3.json', response.json)


def test_find_subject_structural_intersecting_variants_4(client):
    url = find_subject_structural_intersecting_variants_query('subject=HCC1143&ranges=NC_000017.11:43044294-43125364,NC_000013.11:32315507-32400268&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}4.json', response.json)


"""
Find Subject Structural Subsuming Variants Tests
------------------------------------------------
"""

def test_find_subject_structural_subsuming_variants_subject_supplied(client):
    url = find_subject_structural_subsuming_variants_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_structural_subsuming_variants_ranges_supplied(client):
    url = find_subject_structural_subsuming_variants_query('subject=HG00403&testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_structural_subsuming_variants_1(client):
    url = find_subject_structural_subsuming_variants_query('subject=HCC1143&ranges=NC_000007.14:116672195-116798386&specimenIdentifiers=HCC1143-Sp1&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_find_subject_structural_subsuming_variants_2(client):
    url = find_subject_structural_subsuming_variants_query('subject=HCC1143&ranges=NC_000007.14:116672195-116798386&specimenIdentifiers=HCC1143-Sp1&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR}2.json', response.json)


def test_find_subject_structural_subsuming_variants_3(client):
    url = find_subject_structural_subsuming_variants_query('subject=ABC789&ranges=NC_000022.10:42522500-42526812&includeVariants=true')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR}3.json', response.json)


"""
Find Subject Haplotypes Tests
-----------------------------
"""

def test_find_subject_haplotypes_subject_supplied(client):
    url = find_subject_haplotypes_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_haplotypes_genes_supplied(client):
    url = find_subject_haplotypes_query('subject=HG00403&testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_haplotypes_1(client):
    url = find_subject_haplotypes_query('subject=NB6TK328&genes=http://www.genenames.org/geneId|HGNC:2625')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_HAPLOTYPES_OUTPUT_DIR}1.json', response.json)


def test_find_subject_haplotypes_2(client):
    url = find_subject_haplotypes_query('subject=NB6TK328&genes=http://www.genenames.org/geneId|HGNC:2625,http://www.genenames.org/geneId|HGNC:2621')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_HAPLOTYPES_OUTPUT_DIR}2.json', response.json)


def test_find_subject_haplotypes_3(client):
    url = find_subject_haplotypes_query('subject=NB6TK328&genes=http://www.genenames.org/geneId|HGNC:4932')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_HAPLOTYPES_OUTPUT_DIR}3.json', response.json)


def test_find_subject_haplotypes_4(client):
    url = find_subject_haplotypes_query('subject=NB6TK328&genes=HLA-B')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_HAPLOTYPES_OUTPUT_DIR}4.json', response.json)


"""
Find Subject Specific Haplotypes Tests
--------------------------------------
"""

def test_find_subject_specific_haplotypes_subject_supplied(client):
    url = find_subject_specific_haplotypes_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_specific_haplotypes_haplotypes_supplied(client):
    url = find_subject_specific_haplotypes_query('subject=HG00403&testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_specific_haplotypes_1(client):
    url = find_subject_specific_haplotypes_query('subject=XYZ123&haplotypes=CYP2C19')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_SPECIFIC_HAPLOTYPES_OUTPUT_DIR}1.json', response.json)


def test_find_subject_specific_haplotypes_2(client):
    url = find_subject_specific_haplotypes_query('subject=NB6TK328&haplotypes=HLA-B*27')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_SPECIFIC_HAPLOTYPES_OUTPUT_DIR}2.json', response.json)
