import pymongo
import json
import os
import dateutil.parser


def datetime_parser(json_dict):
    for (key, value) in json_dict.items():
        if key=='testDate':
            try:
                json_dict[key] = dateutil.parser.parse(value)
            except:
                pass
    return json_dict

client = pymongo.MongoClient("{MONGO_CLIENT_URL}")

FHIRGenomicsData = client.FHIRGenomicsData

variants = FHIRGenomicsData.Variants

filename = ''
patient_id = ''
with open(filename) as json_file:
    json_object = json.load(json_file, object_hook=datetime_parser)

variants.insert_many(json_object)

with open('patient_tests.json') as js:
    patient_tests = json.load(js, object_hook=datetime_parser)


tests = FHIRGenomicsData.Tests

for test in patient_tests:
    if test['patientID'] == patient_id:
        tests.insert_one(test)
