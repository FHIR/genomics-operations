import tests.utilities as tu


"""
Find Study Metadata Tests
-------------------------
"""


def test_find_study_metadata_subject_supplied(client):
    url = tu.find_study_metadata_query('testDateRange=ge2010-01-01')
    response = client.get(url)

    assert response.status_code == 400


def test_find_study_metadata_1(client):
    url = tu.find_study_metadata_query('subject=ABC123')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_STUDY_METADATA_OUTPUT_DIR}1.json', response.json)


def test_find_study_metadata_2(client):
    url = tu.find_study_metadata_query('subject=HCC1143&ranges=NC_000007.14:55019016-55211628')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_STUDY_METADATA_OUTPUT_DIR}2.json', response.json)


def test_find_study_metadata_3(client):
    url = tu.find_study_metadata_query('subject=HCC1143&ranges=NC_000007.14:55019016-55119016,NC_000007.14:55019016-55211628')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_STUDY_METADATA_OUTPUT_DIR}3.json', response.json)


def test_find_study_metadata_4(client):
    url = tu.find_study_metadata_query('subject=HCC1143&ranges=NC_000006.11:26091303-26091307')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_STUDY_METADATA_OUTPUT_DIR}4.json', response.json)
