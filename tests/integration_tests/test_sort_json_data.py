import json

import tests.utilities as tu

"""
Test JSON Sorting
-------------------------------------
"""
simple = {
    "type": "simple-fruit",
    "example": "apple",
    "seeds?": "yes"
}

simple_nested = {
    "fruits": [
        {
            "type": "simple-fruits",
            "examples": [
                {
                    "name": "apple"
                },
                {
                    "name": "cherries"
                }
            ]
        },
        {
            "type": "aggregate",
            "examples": [
                {
                    "name": "drupelets",
                    "members": [
                        {
                            "common-name": "raspberry"
                        },
                        {
                            "common-name": "blackberry"
                        }
                    ]
                },
                {
                    "name": "achenes",
                    "common-name": "strawberry"
                }
            ]
        }
    ]
}


def test_json_sorting_simple():
    sorted_actual = tu.sort_json_data(simple)
    expected_loc = f'{tu.SORT_JSON_DATA_OUTPUT_DIR}simple.json'
    with open(expected_loc) as expected_file:
        expected_json = json.load(expected_file)
        assert json.dumps(sorted_actual) == json.dumps(expected_json)


def test_json_sorting_nested():
    sorted_actual = tu.sort_json_data(simple_nested)
    expected_loc = f'{tu.SORT_JSON_DATA_OUTPUT_DIR}nested.json'
    with open(expected_loc) as expected_file:
        expected_json = json.load(expected_file)
        assert json.dumps(sorted_actual) == json.dumps(expected_json)
