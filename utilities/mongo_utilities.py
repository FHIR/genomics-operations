import pymongo

# can use this script to help debug mongo queries

# MongoDB Client URIs
FHIR_genomics_data_client_uri = "mongodb+srv://download:download@cluster0.8ianr.mongodb.net/FHIRGenomicsData"
utilities_data_client_uri = "mongodb+srv://download:download@cluster0.8ianr.mongodb.net/UtilitiesData"

# MongoDB Clients
client = pymongo.MongoClient(FHIR_genomics_data_client_uri)
utilities_client = pymongo.MongoClient(utilities_data_client_uri)

# Databases
db = client.FHIRGenomicsData
utilities_db = utilities_client.UtilitiesData


# Collections
variants = db.Variants
molconseqs = db.MolecConseq
tests = db.Tests


# Create records


# READ records
# print(variants.find_one({"patientID":{"$in":["HG00403"]}}))

# print(molconseqs.find_one({"patientID":"TCGA-DD-A1EH","variantID":{"$in": ["0a4f15c4b50b4d3799f223735029bf9c"]}}))

# results=(molconseqs.aggregate([{"$match":{"patientID": "TCGA-DD-A1EH","variantID":{"$in":["0a4f15c4b50b4d3799f223735029bf9c"]}}}]))
# for result in results:
#     print(result)

matchPortion = {'$match':{'patientID':'TCGA-DD-A1EH','SPDI':{'$in':['NC_000001.10:27105549:C:T','NC_000001.10:27106893::T','NC_000007.13:41729606:G:A']}}}
lookupPortion = {'$lookup':{'from':'MolecConseq','let':{'myvariant_id':'$_id'},'pipeline':[{'$match':{'$expr':{'$and':[{'$or':[{'$eq':['$variantID','$$myvariant_id']}]}]}}}],'as':'molecularConsequenceMatches'}}
addFieldsPortion = {'$addFields':{}}
secondMatchPortion = {'$match':{'molecularConsequenceMatches':{'$exists':True,'$not':{'$size':0}}}}
query = [matchPortion,lookupPortion,addFieldsPortion,secondMatchPortion]
print(query)
results=(variants.aggregate(query))
for result in results:
    print(result)


