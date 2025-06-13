import os
import pymongo
from pathlib import Path
from dotenv import load_dotenv
from app import common

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
molCon_db = db.MolecConseq

# Query input parameters
subject = "L2345"
ranges = ["NC_000007.14:55019016-55211628", "NC_000008.11:127735433-127742951", "NC_000002.11:127735431-127742951"]
experimental = True

# Query txImplication collection based on the provided ranges
orGroup = []
for item in ranges:
    refseq, low_high = item.split(':')
    low, high = map(int, low_high.split('-'))
    orElement = {
        '$or': [
            {
                'b38Region.refseq': refseq,
                'b38Region.start': {'$lt': high},
                'b38Region.end': {'$gte': low}
            },
            {
                'b37Region.refseq': refseq,
                'b37Region.start': {'$lt': high},
                'b37Region.end': {'$gte': low}
            }
        ]
    }
    orGroup.append(orElement)
query = {'$or': orGroup}
txImpQueryResults = (txImplication_db.find(query))

# Now to test patient data against the expressions
for txImpResult in txImpQueryResults:
    # First, we gather the patient data
    type = txImpResult["expression"]["constraints"][0]["type"]
    VariantQuery = {
            "$or": [
                {
                    "patientID": subject,
                    "genomicBuild": "GRCh38",
                    "CHROM": txImpResult["b38Region"]["chrom"],
                    "POS": {'$lt': txImpResult["b38Region"]["end"]},
                    "END": {'$gte': txImpResult["b38Region"]["start"]}
                },
                {
                    "patientID": subject,
                    "genomicBuild": "GRCh37",
                    "CHROM": txImpResult["b37Region"]["chrom"],
                    "POS": {'$lt': txImpResult["b37Region"]["end"]},
                    "END": {'$gte': txImpResult["b37Region"]["start"]}
                }
            ]
        }
    variantQueryResults = (variants_db.find(VariantQuery))
    if type in ["DefiningAlleleConstraint"]:
        VariantIDList = []
        for VariantQueryResult in variantQueryResults:
            VariantIDList.append(VariantQueryResult["_id"])
        MolConQuery = {"variantID": {"$in": VariantIDList}}
        molConQueryResults = molCon_db.find(MolConQuery)

    # Now we can compare patient data with txImplication results
    if type in ["CopyCountConstraint"]:
        print("processing CopyCountConstraint...")
        """
        Cycle through each variantQueryResult, checking to see
        if any variant satisfies the expression constraints in txImpResult.
        """
        for variantQueryResult in variantQueryResults:
            criteriaSatisfied = False
            Copies = txImpResult["expression"]["constraints"][0]["copies"]
            if "SPDI" in variantQueryResult:
                REF = (variantQueryResult["SPDI"]).split(":")[2]
                ALT = (variantQueryResult["SPDI"]).split(":")[3]
                if len(REF) > len(ALT) and Copies[1] < 2:
                    criteriaSatisfied = True
            if "CN" in variantQueryResult:
                CN = variantQueryResult["CN"]
                SVTYPE = variantQueryResult["SVTYPE"]
                if CN in Copies:
                    criteriaSatisfied = True
                elif SVTYPE == "DEL" and (0 in Copies or 1 in Copies):
                    criteriaSatisfied = True
            if criteriaSatisfied:
                print(f" - Variant {variantQueryResult['_id']} satisfies txImp {txImpResult['_id']}")
    elif type in ["CopyChangeConstraint"]:
        print("processing CopyChangeConstraint...")
        """
        Cycle through each variantQueryResult, checking to see
        if any variant satisfies the expression constraints in txImpResult.
        """
        for variantQueryResult in variantQueryResults:
            criteriaSatisfied = False
            CopyChange = txImpResult["expression"]["constraints"][0]["copyChange"]
            if "SPDI" in variantQueryResult:
                REF = (variantQueryResult["SPDI"]).split(":")[2]
                ALT = (variantQueryResult["SPDI"]).split(":")[3]
                if len(REF) > len(ALT) and CopyChange == "loss":
                    criteriaSatisfied = True
                    break
            if "CN" in variantQueryResult:
                CN = variantQueryResult["CN"]
                SVTYPE = variantQueryResult["SVTYPE"]
                if (
                    (CN > 2 and CopyChange == "gain")
                    or (CN < 2 and CopyChange == "loss")
                    or (SVTYPE == "DUP" and CopyChange == "gain")
                    or (SVTYPE == "DEL" and CopyChange == "loss")
                ):
                    criteriaSatisfied = True
                    break
            if criteriaSatisfied:
                print(f" - Variant {variantQueryResult['_id']} satisfies txImp {txImpResult['_id']}")
    elif type in ["DefiningAlleleConstraint"]:
        print("processing DefiningAlleleConstraint...")
        """
        Cycle through each molConQueryResult, checking to see
        if any molCon satisfies the expression constraints in txImpResult.
        """
        criteria = txImpResult["expression"]["constraints"][0]["allele"]["expressions"][0]["value"]
        for molConQueryResult in molConQueryResults:
            criteriaSatisfied = False
            variantID = molConQueryResult["variantID"]
            pHGVS = molConQueryResult["pHGVS"]
            if pHGVS == criteria:
                criteriaSatisfied = True
                break
        if criteriaSatisfied:
            print(f" - MolCon {molConQueryResult['_id']} of variant {variantID} satisfies txImp {txImpResult['_id']}")
