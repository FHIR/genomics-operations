import json
import os

from deepdiff import DeepDiff

FIND_SUBJECT_VARIANTS_URL = "/subject-operations/genotype-operations/$find-subject-variants"
FIND_SUBJECT_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_subject_variants/"

FIND_SUBJECT_SPECIFIC_VARIANTS_URL = "/subject-operations/genotype-operations/$find-subject-specific-variants"
FIND_SUBJECT_SPECIFIC_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_subject_specific_variants/"

FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_URL = "/subject-operations/genotype-operations/$find-subject-structural-intersecting-variants"
FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_subject_structural_intersecting_variants/"

FIND_SUBJECT_STRUCTURAL_SUBSUMING_VARIANTS_URL = "/subject-operations/genotype-operations/$find-subject-structural-subsuming-variants"
FIND_SUBJECT_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_subject_structural_subsuming_variants/"

FIND_SUBJECT_HAPLOTYPES_URL = "/subject-operations/genotype-operations/$find-subject-haplotypes"
FIND_SUBJECT_HAPLOTYPES_OUTPUT_DIR = "tests/expected_outputs/find_subject_haplotypes/"

FIND_SUBJECT_SPECIFIC_HAPLOTYPES_URL = "/subject-operations/genotype-operations/$find-subject-specific-haplotypes"
FIND_SUBJECT_SPECIFIC_HAPLOTYPES_OUTPUT_DIR = "tests/expected_outputs/find_subject_specific_haplotypes/"

FIND_SUBJECT_TX_IMPLICATIONS_URL = "/subject-operations/phenotype-operations/$find-subject-tx-implications"
FIND_SUBJECT_TX_IMPLICATIONS_OUTPUT_DIR = "tests/expected_outputs/find_subject_tx_implications/"

FIND_SUBJECT_DX_IMPLICATIONS_URL = "/subject-operations/phenotype-operations/$find-subject-dx-implications"
FIND_SUBJECT_DX_IMPLICATIONS_OUTPUT_DIR = "tests/expected_outputs/find_subject_dx_implications/"

FIND_SUBJECT_MOLEC_CONSEQ_URL = "/subject-operations/phenotype-operations/$find-subject-molecular-consequences"
FIND_SUBJECT_MOLEC_CONSEQ_OUTPUT_DIR = "tests/expected_outputs/find_subject_molecular_consequences/"

FIND_STUDY_METADATA_URL = "/subject-operations/metadata-operations/$find-study-metadata"
FIND_STUDY_METADATA_OUTPUT_DIR = "tests/expected_outputs/find_study_metadata/"

FIND_POPULATION_SPECIFIC_VARIANTS_URL = "/population-operations/genotype-operations/$find-population-specific-variants"
FIND_POPULATION_SPECIFIC_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_population_specific_variants/"

FIND_POPULATION_STRUCTURAL_INTERSECTING_VARIANTS_URL = "/population-operations/genotype-operations/$find-population-structural-intersecting-variants"
FIND_POPULATION_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_population_structural_intersecting_variants/"

FIND_POPULATION_STRUCTURAL_SUBSUMING_VARIANTS_URL = "/population-operations/genotype-operations/$find-population-structural-subsuming-variants"
FIND_POPULATION_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR = "tests/expected_outputs/find_population_structural_subsuming_variants/"

FIND_POPULATION_SPECIFIC_HAPLOTYPES_URL = "/population-operations/genotype-operations/$find-population-specific-haplotypes"
FIND_POPULATION_SPECIFIC_HAPLOTYPES_OUTPUT_DIR = "tests/expected_outputs/find_population_specific_haplotypes/"

FIND_POPULATION_TX_IMPLICATIONS_URL = "/population-operations/phenotype-operations/$find-population-tx-implications"
FIND_POPULATION_TX_IMPLICATIONS_OUTPUT_DIR = "tests/expected_outputs/find_population_tx_implications/"

FIND_POPULATION_DX_IMPLICATIONS_URL = "/population-operations/phenotype-operations/$find-population-dx-implications"
FIND_POPULATION_DX_IMPLICATIONS_OUTPUT_DIR = "tests/expected_outputs/find_population_dx_implications/"

FIND_POPULATION_MOLECULAR_CONSEQUENCES_URL = "/population-operations/phenotype-operations/$find-population-molecular-consequences"
FIND_POPULATION_MOLECULAR_CONSEQUENCES_OUTPUT_DIR = "tests/expected_outputs/find_population_molecular_consequences/"

GET_FEATURE_COORDINATES_URL = "/utilities/get-feature-coordinates"
GET_FEATURE_COORDINATES_OUTPUT_DIR = "tests/expected_outputs/get_feature_coordinates/"

FIND_THE_GENE_URL = "/utilities/find-the-gene"
FIND_THE_GENE_OUTPUT_DIR = "tests/expected_outputs/find_the_gene/"

NORMALIZE_HLA_URL = "/utilities/normalize-hla"
NORMALIZE_HLA_OUTPUT_DIR = "tests/expected_outputs/normalize_hla/"

SORT_JSON_DATA_OUTPUT_DIR = "tests/expected_outputs/sort_json_data/"


def find_subject_variants_query(query):
    return f"{FIND_SUBJECT_VARIANTS_URL}?{query}"


def find_subject_specific_variants_query(query):
    return f"{FIND_SUBJECT_SPECIFIC_VARIANTS_URL}?{query}"


def find_subject_structural_intersecting_variants_query(query):
    return f"{FIND_SUBJECT_STRUCTURAL_INTERSECTING_VARIANTS_URL}?{query}"


def find_subject_structural_subsuming_variants_query(query):
    return f"{FIND_SUBJECT_STRUCTURAL_SUBSUMING_VARIANTS_URL}?{query}"


def find_subject_haplotypes_query(query):
    return f"{FIND_SUBJECT_HAPLOTYPES_URL}?{query}"


def find_subject_specific_haplotypes_query(query):
    return f"{FIND_SUBJECT_SPECIFIC_HAPLOTYPES_URL}?{query}"


def find_subject_tx_implications_query(query):
    return f"{FIND_SUBJECT_TX_IMPLICATIONS_URL}?{query}"


def find_subject_dx_implications_query(query):
    return f"{FIND_SUBJECT_DX_IMPLICATIONS_URL}?{query}"


def find_subject_molecular_consequences_query(query):
    return f"{FIND_SUBJECT_MOLEC_CONSEQ_URL}?{query}"


def find_study_metadata_query(query):
    return f"{FIND_STUDY_METADATA_URL}?{query}"


def find_population_specific_variants_query(query):
    return f"{FIND_POPULATION_SPECIFIC_VARIANTS_URL}?{query}"


def find_population_structural_intersecting_variants_query(query):
    return f"{FIND_POPULATION_STRUCTURAL_INTERSECTING_VARIANTS_URL}?{query}"


def find_population_structural_subsuming_variants_query(query):
    return f"{FIND_POPULATION_STRUCTURAL_SUBSUMING_VARIANTS_URL}?{query}"


def find_population_specific_haplotypes_query(query):
    return f"{FIND_POPULATION_SPECIFIC_HAPLOTYPES_URL}?{query}"


def find_population_tx_implications_query(query):
    return f"{FIND_POPULATION_TX_IMPLICATIONS_URL}?{query}"


def find_population_dx_implications_query(query):
    return f"{FIND_POPULATION_DX_IMPLICATIONS_URL}?{query}"


def find_population_molecular_consequences_query(query):
    return f"{FIND_POPULATION_MOLECULAR_CONSEQUENCES_URL}?{query}"


def get_feature_coordinates_query(query):
    return f"{GET_FEATURE_COORDINATES_URL}?{query}"


def find_the_gene_query(query):
    return f"{FIND_THE_GENE_URL}?{query}"


def normalize_hla_utility_query(query):
    return f"{NORMALIZE_HLA_URL}?{query}"


def sort_json_data(data):
    """
    Recursively sorts the keys of a JSON object or the elements of a JSON list.
    """
    if isinstance(data, dict):
        return {k: sort_json_data(v) for k, v in sorted(data.items())}
    elif isinstance(data, list):
        return sorted((sort_json_data(v) for v in data), key=lambda x: json.dumps(x, sort_keys=True))

    return data


def compare_actual_and_expected_output(filename, actual_json):
    with open(filename) as expected_output_file:
        expected_json = json.load(expected_output_file)
        sorted_expected = sort_json_data(expected_json)
        sorted_actual = sort_json_data(actual_json)

        # TODO: It should suffice to just compare the sorted JSON strings, but this is more flexible
        diff = DeepDiff(sorted_actual, sorted_expected, ignore_order=False, report_repetition=True)
        if diff != {}:
            if 'OVERWRITE_TEST_EXPECTED_DATA' in os.environ:
                with open(filename, 'w') as expected_output_file:
                    json.dump(sorted_actual, expected_output_file, indent=4)
            else:
                assert json.dumps(sorted_expected) == json.dumps(sorted_expected)
