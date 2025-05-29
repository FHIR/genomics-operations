import os
import pymongo
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent.parent / ".env")
# Load secrets from secrets.env file if available
load_dotenv(Path(__file__).parent.parent / "secrets.env")

# MongoDB Client URIs
FHIR_genomics_data_client_uri = f"mongodb+srv://readonly:{os.getenv('MONGODB_READONLY_PASSWORD')}@cluster0.8ianr.mongodb.net/FHIRGenomicsData"

# MongoDB Clients
client = pymongo.MongoClient(FHIR_genomics_data_client_uri)

# Databases
db = client.FHIRGenomicsData

# Collections
patients_db = db.Patients
variants_db = db.Variants
tests_db = db.Tests
genotypes_db = db.Genotypes
dxImplication_db = db.dxImplication
txImplication_db = db.txImplication

# QUERY INPUT PARAMETERS
subject = "L2345"
ranges = ["NC_000007.14:55019016-55211628", "NC_000008.11:127735433-127742951", "NC_000002.11:127735431-127742951"]
experimental = True

refseq = []
low = []
high = []

for item in ranges:
    xrefseq, low_high = item.split(':')
    xlow, xhigh = map(int, low_high.split('-'))
    refseq.append(xrefseq)
    low.append(xlow)
    high.append(xhigh)

orGroup = []
for i in range(len(refseq)):
    item = {
            'region': {
                '$elemMatch': {
                    'refseq': refseq[i],
                    'start': {'$lt': high[i]},
                    'end': {'$gte': low[i]}
                }
            }
        }
    orGroup.append(item)

query = {
    '$or': orGroup
}

resultSet = (txImplication_db.find(query))
for result in resultSet:
    print(result)
