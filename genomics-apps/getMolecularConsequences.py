import streamlit as st
import pandas as pd
import requests
from st_aggrid import GridOptionsBuilder, AgGrid
import pyliftover

st.set_page_config(
    page_title="Get Molecular Consequences",
    page_icon="random",
    layout="wide",
    initial_sidebar_state="collapsed"
)


@st.cache_data
def getFeatureCoordinates(gene):
    url = 'https://fhir-gen-ops.herokuapp.com/utilities/get-feature-coordinates?gene='+gene
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


@st.cache_data
def findSubjectVariants(subject, range, addAnnotationsFlag):
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject=' + \
        subject+'&ranges='+range+'&includeVariants=true'
    if addAnnotationsFlag:
        url += '&includePhasing=true'
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


@st.cache_data
def findSubjectStructuralIntersectingVariants(subject, range, addAnnotationsFlag):
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-structural-intersecting-variants?subject=' + \
        subject+'&ranges='+range+'&includeVariants=true'
    if addAnnotationsFlag:
        url += '&includePhasing=true'
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


@st.cache_data
def findSubjectHaplotypes(subject, geneId):
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-haplotypes?subject='+subject+'&genes='+geneId
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


def _liftOver(chrom, pos):
    lo_b37_to_b38 = pyliftover.LiftOver('hg19', 'hg38')
    newPos = lo_b37_to_b38.convert_coordinate(chrom, pos, '+')
    return {"chrom": chrom, "pos": newPos[0][1]}


@st.cache_data
def findSPDI(SPDI):
    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/spdi/' + \
        SPDI + '/canonical_representative'
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


def computeAnnotations(SPDI):
    # st.write("SPDI:" + SPDI)
    # st.write(findSPDI(SPDI))

    SPDIList = SPDI.split(':')
    chromosome = SPDIList[0].split('.')[0][-2:]
    version = SPDIList[0].split('.')[1]
    # version = int(version) + 1
    version = int(version) + 1

    if chromosome[0] == '0':
        chromosome = chromosome[1]
    chromosome = "chr" + chromosome

    liftOver = _liftOver(chromosome, int(SPDIList[1]))
    # st.write(str(liftOver["pos"]))

    # SPDI = SPDIList[0].split('.')[0] + '.' + str(version) + ":" + str(43049165) + ":" + SPDIList[2] + ":" + SPDIList[3]
    SPDI = SPDIList[0].split('.')[0] + '.' + str(version) + ":" + \
        str(liftOver["pos"]) + ":" + SPDIList[2] + ":" + SPDIList[3]

    # st.write("In computeAnnotations")
    # st.write(SPDI)
    cleanSPDI = SPDI.replace(":", "%3A")

    # st.write(b38SPDI.json())

    # SPDI = SPDIList[0].split('.')[0] + '.' + str(version) + ":" + str(b38SPDI.json()["data"]["position"]) + ":" + SPDIList[2] + ":" + SPDIList[3]

    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/spdi/' + cleanSPDI + '/hgvs'
    headers = {'Accept': 'application/json'}
    HGVS = requests.get(url, headers=headers)

    # st.write(HGVS)
    # st.write(HGVS.json())

    HGVS = HGVS.json()['data']['hgvs']

    server = 'https://rest.ensembl.org'
    ext = '/vep/human/hgvs/' + HGVS
    headers = {"Content-Type": "application/json",
               "Accept": "application/json"}
    r = requests.get(server+ext, headers=headers)

    molecImpact = ""

    # st.write(r.json())

    for entry in r.json()[0]['transcript_consequences']:
        if "impact" in entry:
            molecImpact = "Computed: " + entry["impact"]
            break

    return molecImpact


def findMNVs(cisVariantsID, response, cisVariantsList, convertedVariants):
    for IDList in cisVariantsID:
        # st.write(IDList)
        SPDIList = []
        for i in response.json()["parameter"]:
            for j in i["part"]:
                if j["name"] == "variant":
                    if j["resource"]["id"] in IDList:
                        for k in j["resource"]["component"]:
                            if k["code"]["coding"][0]["code"] == "81252-9" and k["valueCodeableConcept"]["coding"][0]["display"] not in SPDIList:
                                SPDIString = (
                                    k["valueCodeableConcept"]["coding"][0]["display"])
                                SPDIStringList = SPDIString.split(":")
                                SPDIList.append({
                                    "chromosome": SPDIStringList[0],
                                    "position": SPDIStringList[1],
                                    "refAllele": SPDIStringList[2],
                                    "altAllele": SPDIStringList[3]
                                })

                                convertedVariants.append(SPDIString)

        SPDIList = sorted(SPDIList, key=lambda d: d['position'])
        refAllele = ""
        altAllele = ""
        for SPDIDict in SPDIList:
            refAllele += SPDIDict["refAllele"]
            altAllele += SPDIDict["altAllele"]

        # st.write(SPDIList)
        SPDI = SPDIList[0]["chromosome"] + ":" + \
            SPDIList[0]["position"] + ":" + refAllele + ":" + altAllele
        cisVariantsList.append({
            "SPDI": SPDI,
            "Molecular Impact": computeAnnotations(SPDI)
        })


def getImplication(entry):
    for component in entry['resource']['component']:
        if component['code']['coding'][0]['code'] == "53037-8":
            if 'coding' in component['valueCodeableConcept']:
                return component['valueCodeableConcept']['coding'][0]['display']

            if 'text' in component['valueCodeableConcept']:
                return component["valueCodeableConcept"]['text']
            else:
                return "none"

# def putVariant(variantList, entry, pathogenicity):
# 	SPDI = ''
# 	for component in entry['resource']['component']:
# 		if component['code']['coding'][0]['code'] == "81252-9":
# 			SPDI=k["valueCodeableConcept"]["coding"][0]["display"]

# 	variant["Pathogenicity"] = pathogenicity


def findPathogenicities(variantList, subject, range):
    for variant in variantList:

        url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/phenotype-operations/$find-subject-dx-implications?subject=' + \
            subject + '&variants=' + variant['SPDI'].replace(":", "%3A")
        headers = {'Accept': 'application/json'}
        response = requests.get(url, headers=headers)

        if response.status_code != 200:
            # st.write(range)
            continue

        if 'parameter' not in response.json():
            continue

        for parameter in response.json()['parameter']:
            for entry in parameter["part"]:
                if entry["name"] == "implication":
                    variant["Pathogenicity"] = getImplication(entry)

                # if entry["name"] == "variant":
                # 	putVariant(variantList, entry, pathogenicity)


st.title("Get Molecular Consequences")
st.markdown("This app illustrates [FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html) **find-subject-variants**, **find-subject-intersecting-variants**, \
and **find-subject-haplotypes**; and the **get-feature-coordinates** utility. Enter patient and gene in the sidebar and click 'run'. All overlapping simple variants, structural variants, \
and genotypes are returned. Check the compute additional annotations button to calculate annotations for variants that were previously unannotated, concatenate SNVs that are \
in cis into MNVs, and annotate those MNVs.")

with st.sidebar:
    subject = st.text_input("Enter patient ID", value="HG00403")
    gene = st.text_input(
        "Enter gene (HGNC gene symbol or code)", value="BRCA1")
    computeAnnotationsFlag = st.checkbox("Compute additional annotations")

convertedVariants = []
if st.sidebar.button("Run"):
    variantList = []
    cisVariantsID = []
    cisVariantsList = []
    range = {}
    geneFeatures = getFeatureCoordinates(gene)
    if len(geneFeatures.json()) > 0:
        geneId = geneFeatures.json()[0]["geneId"]
        range = geneFeatures.json()[0]["build38Coordinates"]
        # get and parse simple variants
        response = findSubjectVariants(subject, range, computeAnnotationsFlag)
        if response.status_code == 200:
            for i in response.json()["parameter"]:
                for j in i["part"]:
                    if j["name"] == "variant":
                        dnaChangeType = "simple"
                        sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype, molecImpact, copyNumber = "", "", "", "", "", "", "", "", "", "", ""
                        for k in j["resource"]["component"]:
                            if k["code"]["coding"][0]["code"] == "48002-0":
                                sourceClass = k["valueCodeableConcept"]["coding"][0]["display"]
                            elif k["code"]["coding"][0]["code"] == "53034-5":
                                allelicState = k["valueCodeableConcept"]["coding"][0]["display"]
                            elif k["code"]["coding"][0]["code"] == "81252-9":
                                SPDI = k["valueCodeableConcept"]["coding"][0]["display"]
                            elif k["code"]["coding"][0]["code"] == "81258-6":
                                alleleFreq = k["valueQuantity"]["value"]
                            elif k["code"]["coding"][0]["code"] == "molecular-consequence":
                                molecImpact = k["interpretation"]["text"]
                        if molecImpact == "":
                            if computeAnnotationsFlag:
                                molecImpact = computeAnnotations(SPDI)
                        variantList.append({
                            "DNA Change Type": dnaChangeType,
                            "Source Class": sourceClass,
                            "SPDI": SPDI,
                            # "Location": location,
                            # "Genotype": genotype,
                            "Allelic State": allelicState,
                            "Molecular Impact": molecImpact,
                            # "Copy Number": copyNumber,
                            "Allele Frequeny": alleleFreq, })

                    if j["name"] == "sequencePhaseRelationship" and computeAnnotationsFlag:
                        if j["resource"]["valueCodeableConcept"]["coding"][0]["code"] == "Cis":
                            derivedFromIDs = []
                            for idEntry in j["resource"]["derivedFrom"]:
                                derivedFromIDs.append(
                                    idEntry["reference"].split('/')[1])
                            cisVariantsID.append(derivedFromIDs)
            if computeAnnotationsFlag:
                findMNVs(cisVariantsID, response,
                         cisVariantsList, convertedVariants)
        # get and parse structural variants
        response = findSubjectStructuralIntersectingVariants(
            subject, range, computeAnnotationsFlag)
        cisVariantsID = []
        if response.status_code == 200:
            for i in response.json()["parameter"]:
                for j in i["part"]:
                    if j["name"] == "variant":
                        sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype, dnaChangeType, molecImpact, copyNumber = "", "", "", "", "", "", "", "", "", "", "", ""
                        for k in j["resource"]["component"]:
                            if k["code"]["coding"][0]["code"] == "48002-0":
                                sourceClass = k["valueCodeableConcept"]["coding"][0]["display"]
                            elif k["code"]["coding"][0]["code"] == "53034-5":
                                allelicState = k["valueCodeableConcept"]["coding"][0]["display"]
                            elif k["code"]["coding"][0]["code"] == "81258-6":
                                alleleFreq = k["valueQuantity"]["value"]
                            elif k["code"]["coding"][0]["code"] == "48013-7":
                                refSeq = k["valueCodeableConcept"]["coding"][0]["code"]
                            elif k["code"]["coding"][0]["code"] == "82155-3":
                                copyNumber = k["valueQuantity"]["value"]
                            elif k["code"]["coding"][0]["code"] == "48019-4":
                                dnaChangeType = k["valueCodeableConcept"]["coding"][0]["display"]
                            elif k["code"]["coding"][0]["code"] == "81302-2":
                                start = k["valueRange"]["low"]["value"]
                                end = k["valueRange"]["high"]["value"]
                            elif k["code"]["coding"][0]["code"] == "molecular-consequence":
                                molecImpact = k["interpretation"][0]["text"]
                            location = f"{refSeq}:{start}-{end}"
                        if molecImpact == "":
                            "in computemolecImpact"
                            if computeAnnotationsFlag:
                                molecImpact = computeAnnotations(SPDI)

                        variantList.append({
                            "DNA Change Type": dnaChangeType,
                            "Source Class": sourceClass,
                            "SPDI": SPDI,
                            # "Location": location,
                            # "Genotype": genotype,
                            "Allelic State": allelicState,
                            "Molecular Impact": molecImpact,
                            # "Copy Number": copyNumber,
                            "Allele Frequeny": alleleFreq})

                    if j["name"] == "sequencePhaseRelationship" and computeAnnotationsFlag:
                        if j["resource"]["valueCodeableConcept"]["coding"][0]["code"] == "Cis":
                            derivedFromIDs = []
                            for idEntry in j["resource"]["derivedFrom"]:
                                derivedFromIDs.append(idEntry["reference"])
                            cisVariantsID.append(derivedFromIDs)
            if computeAnnotationsFlag:
                findMNVs(cisVariantsID, response,
                         cisVariantsList, convertedVariants)
        # get and parse haplotypes
        response = findSubjectHaplotypes(subject, geneId)
        if response.status_code == 200:
            for i in response.json()["parameter"]:
                for j in i["part"]:
                    if j["name"] == "genotype":
                        dnaChangeType = "genotype"
                        sourceClass, allelicState, SPDI, alleleFreq, location, refSeq, start, end, genotype, copyNumber = "", "", "", "", "", "", "", "", "", ""
                        genotype = j["resource"]["valueCodeableConcept"]["coding"][0]["display"]
                        variantList.append({
                            "DNA Change Type": dnaChangeType,
                            "Source Class": sourceClass,
                            "SPDI": SPDI,
                            # "Location": location,
                            # "Genotype": genotype,
                            "Allelic State": allelicState,
                            # "Copy Number": copyNumber,
                            "Allele Frequeny": alleleFreq})

    # findPathogenicities(variantList, subject, range)

    data = (pd.DataFrame(variantList))

    gb = GridOptionsBuilder.from_dataframe(data)
    highlightedRows = []
    for index, row in data.iterrows():
        if row['SPDI'] in convertedVariants:
            highlightedRows.append(index)

    gb.configure_selection('multiple', pre_selected_rows=highlightedRows)
    AgGrid(data, enable_enterprise_modules=True, update_mode="value_changed",
           allow_unsafe_jscode=True, gridOptions=gb.build())

    # data.apply(colors.get, subset=['cisVariantIndex'])

    # AgGrid(data.style.format('{:.0f}').applymap(colors.get, subset=['cisVariantIndex']).hide_columns(['cisVariantIndex']), \

    # AgGrid(data, \
    # 	enable_enterprise_modules=True, update_mode="value_changed", allow_unsafe_jscode=True)
    st.download_button("Download table (json)",
                       data.T.to_json(), mime="application/json")
    st.download_button("Download table (csv)", data.to_csv(), mime="text/csv")

    findPathogenicities(cisVariantsList, subject, range)

    cisVariantsTable = pd.DataFrame(cisVariantsList)
    mnvgb = GridOptionsBuilder.from_dataframe(cisVariantsTable)
    # AgGrid(cisVariantsTable.style.format('{:.0f}').apply(colors.get, subset=['cisVariantIndex']).hide('cisVariantIndex'), \
    # AgGrid(cisVariantsTable.style.format('{:.0f}').apply(colors.get, subset=['cisVariantIndex']), \
    # 	enable_enterprise_modules=True, update_mode="value_changed", allow_unsafe_jscode=True, key="cisVariants")
    if computeAnnotationsFlag:
        st.write('Concatenated MNV Table')

        mnvgb.configure_selection('multiple', pre_selected_rows=['0'])
        AgGrid(cisVariantsTable, enable_enterprise_modules=True,
               update_mode="value_changed", allow_unsafe_jscode=True, gridOptions=mnvgb.build())

        st.download_button("Download table (json)",
                           cisVariantsTable.T.to_json(), mime="application/json")
        st.download_button("Download table (csv)",
                           cisVariantsTable.to_csv(), mime="text/csv")
