import tests.utilities as tu


"""
Find The Gene Utility Tests
---------------------------
"""


def test_get_find_the_gene_parameter_supplied(client):
    url = tu.find_the_gene_query('')
    response = client.get(url)

    assert response.status_code == 400


def test_get_find_the_gene_1(client):
    url = tu.find_the_gene_query('range=NC_000001.11:11794399-11794400')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_THE_GENE_OUTPUT_DIR}1.json', response.json)
