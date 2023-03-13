import tests.utilities as tu


"""
Find Population Tx Implications Tests
-------------------------------------
"""

# def test_find_population_tx_implications_1(client):
#     url = find_population_tx_implications_query('variants=NC_000006.11:26091178:C:G&includePatientList=true')
#     response = client.get(url)

#     compare_actual_and_expected_output(f'{FIND_POPULATION_TX_IMPLICATIONS_OUTPUT_DIR}1.json', response.json)


def test_find_population_tx_implications_2(client):
    url = tu.find_population_tx_implications_query('treatments=http://www.nlm.nih.gov/research/umls/rxnorm|704&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_TX_IMPLICATIONS_OUTPUT_DIR}2.json', response.json)


def test_find_population_tx_implications_3(client):
    url = tu.find_population_tx_implications_query('treatments=http://www.nlm.nih.gov/research/umls/rxnorm|1147220&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_TX_IMPLICATIONS_OUTPUT_DIR}3.json', response.json)


def test_find_population_tx_implications_4(client):
    url = tu.find_population_tx_implications_query('haplotypes=CYP2C19 *1/*17,CYP2D6 *1/*41&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_TX_IMPLICATIONS_OUTPUT_DIR}4.json', response.json)


def test_find_population_tx_implications_5(client):
    url = tu.find_population_tx_implications_query('conditions=https://disease-ontology.org|3908&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_TX_IMPLICATIONS_OUTPUT_DIR}5.json', response.json)


"""
Find Population Dx Implications Tests
-------------------------------------
"""


def test_find_population_dx_implications_1(client):
    url = tu.find_population_dx_implications_query('conditions=https://www.ncbi.nlm.nih.gov/medgen|C3469186&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_DX_IMPLICATIONS_OUTPUT_DIR}1.json', response.json)


def test_find_population_dx_implications_2(client):
    url = tu.find_population_dx_implications_query('conditions=https://www.ncbi.nlm.nih.gov/medgen|C1708353')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_DX_IMPLICATIONS_OUTPUT_DIR}2.json', response.json)


def test_find_population_dx_implications_3(client):
    url = tu.find_population_dx_implications_query('conditions=https://www.ncbi.nlm.nih.gov/medgen|C3469186&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_DX_IMPLICATIONS_OUTPUT_DIR}3.json', response.json)


def test_find_population_dx_implications_4(client):
    url = tu.find_population_dx_implications_query(
        'conditions=https://www.ncbi.nlm.nih.gov/medgen|C0677776,https://www.ncbi.nlm.nih.gov/medgen|C4552100,https://www.ncbi.nlm.nih.gov/medgen|C0020445&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_DX_IMPLICATIONS_OUTPUT_DIR}4.json', response.json)
