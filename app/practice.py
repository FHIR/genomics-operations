import os
import pymongo
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent.parent / ".env")
load_dotenv(Path(__file__).parent.parent / "secrets.env")

FHIR_genomics_data_client_uri = f"mongodb+srv://readonly:{os.getenv('MONGODB_READONLY_PASSWORD')}@cluster0.8ianr.mongodb.net/FHIRGenomicsData"
client = pymongo.MongoClient(FHIR_genomics_data_client_uri)
db = client.FHIRGenomicsData

patients_db = db.Patients
variants_db = db.Variants
tests_db = db.Tests
genotypes_db = db.Genotypes
dxImplication_db = db.dxImplication
txImplication_db = db.txImplication

# ...existing code...

subject = "L2345"
ranges = [
    "NC_000007.14:55019016-55211628",
    "NC_000008.11:127735433-127742951",
    "NC_000002.11:127735431-127742951"
]
experimental = True

orGroup = []
for item in ranges:
    refseq, low_high = item.split(':')
    low, high = map(int, low_high.split('-'))
    orGroup.append({
        'b38Region.refseq': refseq,
        'b38Region.start': {'$lt': high},
        'b38Region.end': {'$gte': low}
    })
    orGroup.append({
        'b37Region.refseq': refseq,
        'b37Region.start': {'$lt': high},
        'b37Region.end': {'$gte': low}
    })

query = {'$or': orGroup}

resultSet = txImplication_db.find(query)
# ...existing code...

unique_lines = set()

for result in resultSet:
    constraint_type = result["expression"]["constraints"][0]["type"]
    if constraint_type == "CopyChangeConstraint":
        refseq = result["expression"]["constraints"][1]["location"]["sequenceReference"]["id"].split(":")[1]
        start = result["expression"]["constraints"][1]["location"]["start"]
        end = result["expression"]["constraints"][1]["location"]["end"]
    elif constraint_type == "CopyCountConstraint":
        refseq = result["expression"]["constraints"][1]["location"]["sequenceReference"]["id"].split(":")[1]
        start = result["expression"]["constraints"][1]["location"]["start"]
        end = result["expression"]["constraints"][1]["location"]["end"]
    elif constraint_type == "DefiningAlleleConstraint":
        refseq = result["expression"]["constraints"][0]["allele"]["location"]["sequenceReference"]["id"].split(":")[1]
        start = result["expression"]["constraints"][0]["allele"]["location"]["start"]
        end = result["expression"]["constraints"][0]["allele"]["location"]["end"]
    else:
        refseq, start, end = None, None, None

    line = f"Type: {constraint_type}, RefSeq: {refseq}, Start: {start}, End: {end}"
    unique_lines.add(line)

for line in unique_lines:
    print(line)
# ...existing code...
