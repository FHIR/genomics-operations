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

subject = "L2345"
ranges = ["NC_000007.14:55019016-55211628", "NC_000008.11:127735433-127742951"]
experimental = True

# Query mongoDb for txImplication records that have regions that intersect our ranges

# print(patients_db.find_one({"patientID": subject}))
# print(txImplication_db.find_one({"evidenceLevel": "CPIC Level A"}))
resultSet = (txImplication_db.find({"region": "NC_000007.14:55019016-55211628"}))
for result in resultSet:
    print(result)
