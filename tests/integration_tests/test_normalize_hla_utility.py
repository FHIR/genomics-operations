import tests.utilities as tu

"""
Normalize HLA Utility Tests
--------------------------------
"""


def test_normalize_hla_utility(client):
    url = tu.normalize_hla_utility_query('allele=B14')
    response = client.get(url)

    assert response.status_code == 200

    # This utility is not deterministic for a few output items such as `exon` and `U2`
    # tu.compare_actual_and_expected_output(f'{tu.NORMALIZE_HLA_OUTPUT_DIR}1.json', response.json)
