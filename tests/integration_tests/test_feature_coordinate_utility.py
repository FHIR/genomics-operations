import pytest
from tests.utilities import *


"""
Feature Coordinate Utility Tests
--------------------------------
"""

def test_get_feature_coordinates_one_parameter_supplied(client):
    url = get_feature_coordinates_query('')
    response = client.get(url)

    assert response.status_code == 400


def test_get_feature_coordinates_only_one_parameter_supplied(client):
    url = get_feature_coordinates_query('chromosome=chr1,transcript=NM_000014.6')
    response = client.get(url)

    assert response.status_code == 400


def test_get_feature_coordinates_1(client):
    url = get_feature_coordinates_query('chromosome=chr1')
    response = client.get(url)

    compare_actual_and_expected_output(f'{GET_FEATURE_COORDINATES_OUTPUT_DIR}1.json', response.json)


def test_get_feature_coordinates_2(client):
    url = get_feature_coordinates_query('chromosome=chrM')
    response = client.get(url)

    compare_actual_and_expected_output(f'{GET_FEATURE_COORDINATES_OUTPUT_DIR}2.json', response.json)


def test_get_feature_coordinates_3(client):
    url = get_feature_coordinates_query('chromosome=chry')
    response = client.get(url)

    compare_actual_and_expected_output(f'{GET_FEATURE_COORDINATES_OUTPUT_DIR}3.json', response.json)


def test_get_feature_coordinates_4(client):
    url = get_feature_coordinates_query('gene=HGNC:5')
    response = client.get(url)

    compare_actual_and_expected_output(f'{GET_FEATURE_COORDINATES_OUTPUT_DIR}4.json', response.json)


def test_get_feature_coordinates_5(client):
    url = get_feature_coordinates_query('transcript=NM_000014.6')
    response = client.get(url)

    compare_actual_and_expected_output(f'{GET_FEATURE_COORDINATES_OUTPUT_DIR}5.json', response.json)


def test_get_feature_coordinates_6(client):
    url = get_feature_coordinates_query('protein=NP_000005.3')
    response = client.get(url)

    compare_actual_and_expected_output(f'{GET_FEATURE_COORDINATES_OUTPUT_DIR}6.json', response.json)
