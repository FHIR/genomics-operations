import streamlit as st
import pandas as pd
import requests
import json
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode
import string
import pyliftover

st.set_page_config(
     page_title="Get Variants",
     page_icon="random",
     layout="wide",
     initial_sidebar_state="collapsed"
 )

@st.cache
def getFeatureCoordinates(gene):
	url='https://fhir-gen-ops.herokuapp.com/utilities/get-feature-coordinates?gene='+gene
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

@st.cache
def findSubjectVariants(subject, range, addAnnotationsFlag):
	url='https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='+subject+'&ranges='+range+'&includeVariants=true'
	if addAnnotationsFlag:
		url += '&includePhasing=true'
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

@st.cache
def findSubjectStructuralIntersectingVariants(subject, range, addAnnotationsFlag):
	url='https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-structural-intersecting-variants?subject='+subject+'&ranges='+range+'&includeVariants=true'
	if addAnnotationsFlag:
		url += '&includePhasing=true'
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

@st.cache
def findSubjectHaplotypes(subject, geneId):
	url='https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-haplotypes?subject='+subject+'&genes='+geneId
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

def _liftOver(chrom, pos):
	lo_b37_to_b38 = pyliftover.LiftOver('hg19','hg38')
	newPos = lo_b37_to_b38.convert_coordinate(chrom,pos,'+')
	return {"chrom":chrom,"pos":newPos[0][1]}

def computeAnnotations(SPDI):
	st.write("SPDI:" + SPDI)
	SPDIList = SPDI.split(':')
	chromosome = SPDIList[0].split('.')[0][-2]
	version = SPDIList[0].split('.')[1]
	version = int(version) + 1
	if chromosome[0] == '0':
		chromosome = chromosome[1]
	chromosome = "chr" + chromosome
	liftOver = _liftOver(chromosome, int(SPDIList[1]))

	# SPDI = SPDIList[0].split('.')[0] + '.' + str(version) + ":" + str(43049165) + ":" + SPDIList[2] + ":" + SPDIList[3]
	SPDI = SPDIList[0].split('.')[0] + '.' + str(version) + ":" + str(liftOver["pos"]) + ":" + SPDIList[2] + ":" + SPDIList[3]

	st.write("In computeAnnotations")
	st.write(SPDI)
	cleanSPDI = SPDI.replace(":", "%3A")
	st.write(cleanSPDI)
	url='https://api.ncbi.nlm.nih.gov/variation/v0/spdi/' + cleanSPDI + '/hgvs'
	headers={'Accept': 'application/json'}
	HGVS = requests.get(url, headers=headers)

	st.write(HGVS)
	st.write(HGVS.json())

	HGVS = HGVS.json()['data']['hgvs']

	server='https://rest.ensembl.org'
	ext = '/vep/human/hgvs/' + HGVS
	headers={ "Content-Type" : "application/json", "Accept" : "application/json" }
	r = requests.get(server+ext, headers=headers)

	molecImpact = ""

	st.write(r.json())

	for entry in r.json()[0]['transcript_consequences']:
		if "impact" in entry:
			molecImpact = "Computed: " + entry["impact"]
			break

	return molecImpact

def findMNVs(cisVariantsID, response, cisVariantsList):
	for IDList in cisVariantsID:
		st.write(IDList)
		SPDIList = []
		for i in response.json()["parameter"]:
					for j in i["part"]:
						if j["name"]=="variant":
							if j["resource"]["id"] in IDList:
								for k in j["resource"]["component"]:
									if k["code"]["coding"][0]["code"]=="81252-9" and k["valueCodeableConcept"]["coding"][0]["display"] not in SPDIList:
										SPDIString = (k["valueCodeableConcept"]["coding"][0]["display"])
										SPDIString = SPDIString.split(":")
										SPDIList.append({
											"chromosome": SPDIString[0],
											"position": SPDIString[1],
											"refAllele": SPDIString[2],
											"altAllele": SPDIString[3]
										})
				
		SPDIList = sorted(SPDIList, key=lambda d: d['position'])
		refAllele = ""
		altAllele = ""
		for SPDIDict in SPDIList:
			refAllele += SPDIDict["refAllele"]
			altAllele += SPDIDict["altAllele"]
		
		st.write(SPDIList)
		SPDI = SPDIList[0]["chromosome"] + ":" + SPDIList[0]["position"] + ":" + refAllele + ":" + altAllele
		cisVariantsList.append({
			"SPDI": SPDI,
			"Molecular Impact": computeAnnotations(SPDI)
		})
	
	return cisVariantsList



st.title("Get Variants")
st.markdown("This app illustrates [FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html) find-subject-variants, find-subject-intersecting-variants, \
	and find-subject-haplotypes; and the get-feature-coordinates utility. Enter patient and gene in the sidebar and click 'run'. All overlapping simple variants, structural variants, \
	and genotypes are returned.")

with st.sidebar:
	subject=st.text_input("Enter patient ID", value="HG00403")
	gene=st.text_input("Enter gene (HGNC gene symbol or code)",value="APC")
	computeAnnotationsFlag = st.checkbox("Compute additional annotations")

if st.sidebar.button("Run"):
	variantList=[]
	cisVariantsID = []
	cisVariantsList = []
	geneFeatures=getFeatureCoordinates(gene)
	if len(geneFeatures.json())>0:
		geneId=geneFeatures.json()[0]["geneId"]
		range=geneFeatures.json()[0]["build38Coordinates"]
		# get and parse simple variants
		response=findSubjectVariants(subject,range, computeAnnotationsFlag)
		if response.status_code == 200:
			for i in response.json()["parameter"]:
				for j in i["part"]:
					if j["name"]=="variant":
						dnaChangeType="simple"
						sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype , molecImpact, copyNumber="","", "", "", "", "", "", "", "", "", ""
						for k in j["resource"]["component"]:
							if k["code"]["coding"][0]["code"]=="48002-0":
								sourceClass=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="53034-5":
								allelicState=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="81252-9":
								SPDI=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="81258-6":
								alleleFreq=k["valueQuantity"]["value"]
							elif k["code"]["coding"][0]["code"]=="molecular-consequence":
								k["interpretation"]
								molecImpact = k["interpretation"]["text"]
						if molecImpact == "":
							if computeAnnotationsFlag:
								molecImpact = computeAnnotations(SPDI)
						variantList.append({
							"DNA Change Type": dnaChangeType,
							"Source Class": sourceClass,
							"SPDI": SPDI,
							"Location": location,
							"Genotype": genotype,
							"Allelic State": allelicState,
							"Molecular Impact": molecImpact,
							"Copy Number": copyNumber,
							"Allele Frequeny": alleleFreq})

					if j["name"] == "sequencePhaseRelationship" and computeAnnotationsFlag:
						if j["resource"]["valueCodeableConcept"]["coding"][0]["code"] == "Cis":
							derivedFromIDs = []
							for idEntry in j["resource"]["derivedFrom"] :
								derivedFromIDs.append(idEntry["reference"].split('/')[1])
							cisVariantsID.append(derivedFromIDs)
			if computeAnnotationsFlag:
				cisVariantsList = findMNVs(cisVariantsID, response, cisVariantsList)
		# get and parse structural variants
		response=findSubjectStructuralIntersectingVariants(subject,range, computeAnnotationsFlag)
		cisVariantsID = []
		if response.status_code == 200:
			for i in response.json()["parameter"]:
				for j in i["part"]:
					if j["name"]=="variant":
						sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype, dnaChangeType , molecImpact, copyNumber="","", "", "", "", "", "", "", "", "", "", ""
						for k in j["resource"]["component"]:
							if k["code"]["coding"][0]["code"]=="48002-0":
								sourceClass=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="53034-5":
								allelicState=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="81258-6":
								alleleFreq=k["valueQuantity"]["value"]
							elif k["code"]["coding"][0]["code"]=="48013-7":
								refSeq=k["valueCodeableConcept"]["coding"][0]["code"]
							elif k["code"]["coding"][0]["code"]=="82155-3":
								copyNumber=k["valueQuantity"]["value"]
							elif k["code"]["coding"][0]["code"]=="48019-4":
								dnaChangeType=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="81302-2":
								start=k["valueRange"]["low"]["value"]								
								end=k["valueRange"]["high"]["value"]
							elif k["code"]["coding"][0]["code"]=="molecular-consequence":
								molecImpact = k["interpretation"][0]["text"]
							location=f"{refSeq}:{start}-{end}"
						if molecImpact == "":
							"in computemolecImpact"
							if computeAnnotationsFlag:
								molecImpact = computeAnnotations(SPDI)

						variantList.append({
							"DNA Change Type": dnaChangeType,
							"Source Class": sourceClass,
							"SPDI": SPDI,
							"Location": location,
							"Genotype": genotype,
							"Allelic State": allelicState,
							"Molecular Impact": molecImpact,
							"Copy Number": copyNumber,
							"Allele Frequeny": alleleFreq})

					if j["name"] == "sequencePhaseRelationship" and computeAnnotationsFlag:
						if j["resource"]["valueCodeableConcept"]["coding"][0]["code"] == "Cis":
							derivedFromIDs = []
							for idEntry in j["resource"]["derivedFrom"] :
								derivedFromIDs.append(idEntry["reference"])
							cisVariantsID.append(derivedFromIDs)
			if computeAnnotationsFlag:
				cisVariantsList = findMNVs(cisVariantsID, response, cisVariantsList)
		# get and parse haplotypes
		response=findSubjectHaplotypes(subject,geneId)
		if response.status_code == 200:
			for i in response.json()["parameter"]:
				for j in i["part"]:
					if j["name"]=="genotype":
						dnaChangeType="genotype"
						sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype, copyNumber="", "", "", "", "", "", "", "", "", ""
						genotype=j["resource"]["valueCodeableConcept"]["coding"][0]["display"]
						variantList.append({
							"DNA Change Type": dnaChangeType,
							"Source Class": sourceClass,
							"SPDI": SPDI,
							"Location": location,
							"Genotype": genotype,
							"Allelic State": allelicState,
							"Copy Number": copyNumber,
							"Allele Frequeny": alleleFreq})

	data=(pd.DataFrame(variantList))
	AgGrid(data, enable_enterprise_modules=True, update_mode="value_changed", allow_unsafe_jscode=True)
	st.download_button("Download table (json)",data.T.to_json(),mime="application/json")
	st.download_button("Download table (csv)",data.to_csv(),mime="text/csv")

	cisVariantsTable = pd.DataFrame(cisVariantsList)
	AgGrid(cisVariantsTable, enable_enterprise_modules=True, update_mode="value_changed", allow_unsafe_jscode=True, key="cisVariants")