import streamlit as st
import pandas as pd
import requests
import json
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode

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
def findSubjectVariants(subject, range):
	url='https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='+subject+'&ranges='+range+'&includeVariants=true'
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

@st.cache
def findSubjectStructuralIntersectingVariants(subject, range):
	url='https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-structural-intersecting-variants?subject='+subject+'&ranges='+range+'&includeVariants=true'
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

@st.cache
def findSubjectHaplotypes(subject, geneId):
	url='https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-haplotypes?subject='+subject+'&genes='+geneId
	headers={'Accept': 'application/json'}
	return requests.get(url, headers=headers)

st.title("Get Variants")
st.markdown("This app illustrates [FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html) find-subject-variants, find-subject-intersecting-variants, \
	and find-subject-haplotypes; and the get-feature-coordinates utility. Enter patient and gene in the sidebar and click 'run'. All overlapping simple variants, structural variants, \
	and genotypes are returned.")

with st.sidebar:
	subject=st.text_input("Enter patient ID", value="HG00403")
	gene=st.text_input("Enter gene (HGNC gene symbol or code)",value="APC")

if st.sidebar.button("Run"):
	variantList=[]
	geneFeatures=getFeatureCoordinates(gene)
	if len(geneFeatures.json())>0:
		geneId=geneFeatures.json()[0]["geneId"]
		range=geneFeatures.json()[0]["build38Coordinates"]
		# get and parse simple variants
		response=findSubjectVariants(subject,range)
		if response.status_code == 200:
			for i in response.json()["parameter"]:
				for j in i["part"]:
					if j["name"]=="variant":
						dnaChangeType="simple"
						sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype , copyNumber="","", "", "", "", "", "", "", "", ""
						for k in j["resource"]["component"]:
							if k["code"]["coding"][0]["code"]=="48002-0":
								sourceClass=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="53034-5":
								allelicState=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="81252-9":
								SPDI=k["valueCodeableConcept"]["coding"][0]["display"]
							elif k["code"]["coding"][0]["code"]=="81258-6":
								alleleFreq=k["valueQuantity"]["value"]
						variantList.append({
							"DNA Change Type": dnaChangeType,
							"Source Class": sourceClass,
							"SPDI": SPDI,
							"Location": location,
							"Genotype": genotype,
							"Allelic State": allelicState,
							"Copy Number": copyNumber,
							"Allele Frequeny": alleleFreq})
		# get and parse structural variants
		response=findSubjectStructuralIntersectingVariants(subject,range)
		if response.status_code == 200:
			for i in response.json()["parameter"]:
				for j in i["part"]:
					if j["name"]=="variant":
						sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype, dnaChangeType , copyNumber="","", "", "", "", "", "", "", "", "", ""
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
							location=f"{refSeq}:{start}-{end}"
						variantList.append({
							"DNA Change Type": dnaChangeType,
							"Source Class": sourceClass,
							"SPDI": SPDI,
							"Location": location,
							"Genotype": genotype,
							"Allelic State": allelicState,
							"Copy Number": copyNumber,
							"Allele Frequeny": alleleFreq})
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
	AgGrid(data)
	st.download_button("Download table (json)",data.T.to_json(),mime="application/json")
	st.download_button("Download table (csv)",data.to_csv(),mime="text/csv")
