import pytest
from tests.utilities import *


"""
Find Subject Tx Implications Tests
----------------------------------
"""

def test_find_subject_tx_implications_subject_supplied(client):
    url = find_subject_tx_implications_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_tx_implications_1(client):
    url = find_subject_tx_implications_query('conditions=https://disease-ontology.org|3908,https://disease-ontology.org|1324,https://disease-ontology.org|3910&subject=CA12345&variants=NM_001354609.2:c.1799T>A')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_TX_IMPLICATIONS_OUTPUT_DIR}1.json', response.json)


def test_find_subject_tx_implications_2(client):
    url = find_subject_tx_implications_query('subject=NB6TK328&treatments=http://www.nlm.nih.gov/research/umls/rxnorm|704')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_TX_IMPLICATIONS_OUTPUT_DIR}2.json', response.json)


def test_find_subject_tx_implications_3(client):
    url = find_subject_tx_implications_query('subject=CA12345&treatments=http://www.nlm.nih.gov/research/umls/rxnorm|1147220')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_TX_IMPLICATIONS_OUTPUT_DIR}3.json', response.json)


def test_find_subject_tx_implications_4(client):
    url = find_subject_tx_implications_query('subject=NB6TK328&haplotypes=CYP2C19 *1/*17,CYP2D6 *1/*41')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_TX_IMPLICATIONS_OUTPUT_DIR}4.json', response.json)


def test_find_subject_tx_implications_5(client):
    url = find_subject_tx_implications_query('subject=CA12345&variants=NM_002524.5:c.182A>C,NM_001354609.2:c.1799T>A&conditions=https://disease-ontology.org|3908')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_TX_IMPLICATIONS_OUTPUT_DIR}5.json', response.json)


"""
Find Subject Dx Implications Tests
----------------------------------
"""

def test_find_subject_dx_implications_subject_supplied(client):
    url = find_subject_dx_implications_query('testDateRange=ge2010-01-01&includeVariants=true')
    response = client.get(url)

    assert response.status_code == 400


def test_find_subject_dx_implications_1(client):
    url = find_subject_dx_implications_query('subject=huC30902&variants=NC_000001.10:161333381:C:T')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_DX_IMPLICATIONS_OUTPUT_DIR}1.json', response.json)


def test_find_subject_dx_implications_2(client):
    url = find_subject_dx_implications_query('subject=huC30902&conditions=https://www.ncbi.nlm.nih.gov/medgen|C1708353')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_DX_IMPLICATIONS_OUTPUT_DIR}2.json', response.json)


def test_find_subject_dx_implications_3(client):
    url = find_subject_dx_implications_query('subject=NB6TK329&conditions=https://www.ncbi.nlm.nih.gov/medgen|C0677776,https://www.ncbi.nlm.nih.gov/medgen|C4552100,https://www.ncbi.nlm.nih.gov/medgen|C0020445')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_DX_IMPLICATIONS_OUTPUT_DIR}3.json', response.json)


def test_find_subject_dx_implications_4(client):
    url = find_subject_dx_implications_query('subject=HG02657&conditions=https://www.ncbi.nlm.nih.gov/medgen|C3469186')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_DX_IMPLICATIONS_OUTPUT_DIR}4.json', response.json)

def test_find_subject_dx_implications_5(client):
    """
    Query with multiple ranges which have different chromosomes
    """
    url = find_subject_dx_implications_query('subject=HG00403&ranges=NC_000001.11:237042183-237833988,NC_000005.9:112043194-112181936,NC_000019.9:11200138-11244496,NC_000013.10:32889644-32974405,NC_000014.9:23412739-23435677')
    response = client.get(url)

    compare_actual_and_expected_output(f'{FIND_SUBJECT_DX_IMPLICATIONS_OUTPUT_DIR}5.json', response.json)
