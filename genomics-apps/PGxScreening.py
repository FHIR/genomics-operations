import streamlit as st
import pandas as pd
import requests
import csv
import random
from st_aggrid import GridOptionsBuilder, AgGrid, JsCode

st.set_page_config(
    page_title="PGx Screening",
    page_icon="random",
    layout="wide",
    initial_sidebar_state="collapsed"
)


@st.cache
def findSubjectHaplotypes(subject):
    # CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP3A5, NUDT15, SLCO1B1, TPMP, UGT1A1
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-haplotypes?subject='+subject + \
        '&genes=http://www.genenames.org/geneId|HGNC:2615,http://www.genenames.org/geneId|HGNC:2623,http://www.genenames.org/geneId|HGNC:2621,http://www.genenames.org/geneId|HGNC:2625,http://www.genenames.org/geneId|HGNC:2638,http://www.genenames.org/geneId|HGNC:23063,http://www.genenames.org/geneId|HGNC:10959,http://www.genenames.org/geneId|HGNC:12014,http://www.genenames.org/geneId|HGNC:12530'
    headers = {'Accept': 'application/json'}
    r = requests.get(url, headers=headers)
    return r.json()


@st.cache
def findSubjectTxImplications(subject, haplotypes):
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/phenotype-operations/$find-subject-tx-implications?subject='+subject+'&haplotypes='+haplotypes
    headers = {'Accept': 'application/json'}
    r = requests.get(url, headers=headers)
    return r.json()


def getMedicationList(subject):
    with open('genomics-apps/data/product.csv') as productFile:
        products = csv.reader(productFile, delimiter=',', quotechar='"')
        productList = []
        for row in products:
            productList.append({"ing": row[0], "product": row[3]})
    medList = []
    for i in range(1, random.randrange(5, 50)):
        j = random.randrange(1, 1383)
        medList.append(productList[j])
    ingList = []
    finalMedList = []
    for i in medList:
        if i["ing"] not in ingList:
            ingList.append(i["ing"])
            finalMedList.append(i)
    return finalMedList


st.title("PGx Screening")
st.markdown("This app illustrates [FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html) find-subject-haplotypes and \
     find-subject-tx-implications. Select patient from dropdown. List of potential drug-gene interactions are shown in right column. \
     The patient's med list is shown in left column, where those meds with PGx interactions are flagged.")


subject = st.selectbox("Select subject:", ('XYZ234', 'NB6TK328', 'NB6TK329', 'XYZ123', 'XYZ234', 'XYZ345'))

col1, col2 = st.columns([2, 3])

genotypeList = ""
patientList = []
evidenceLevel = []
medicationCode = []
medicationName = []
implication = []
genotype = []
PharmGKBID = []

# call findSubjectHaplotypes
haplotypes = findSubjectHaplotypes(subject)

# extract haplotypes
for i in haplotypes["parameter"]:
    for j in i["part"]:
        if j["name"] == "genotype":
            genotypeList = genotypeList+j["resource"]["valueCodeableConcept"]["coding"][0]["code"]+","

# call findSubjectTxImplications
implications = findSubjectTxImplications(subject, genotypeList[:-1])

# extract implications
try:
    count = len(implications["parameter"])
except KeyError:
    count = 0

if count > 0:
    for i in implications["parameter"]:
        for j in i["part"]:
            if j["name"] == "implication":
                # get evidenceLevel, medicationCode, medicationName, implication
                for l in j["resource"]["component"]:
                    if l["code"]["coding"][0]["code"] == "93044-6":
                        evidenceLevel.append(l["valueCodeableConcept"]["text"])
                    elif l["code"]["coding"][0]["code"] == "51963-7":
                        medicationCode.append(l["valueCodeableConcept"]["coding"][0]["code"])
                        medicationName.append(l["valueCodeableConcept"]["coding"][0]["display"])
                    elif l["code"]["coding"][0]["code"] == "predicted-therapeutic-implication":
                        implication.append(l["valueCodeableConcept"]["text"])
            elif j["name"] == "genotype":
                # get genotype, PharmGKBID
                genotype.append(j["resource"]["valueCodeableConcept"]["coding"][0]["code"])
                for l in j["resource"]["component"]:
                    if l["code"]["coding"][0]["code"] == "81252-9":
                        PharmGKBID.append("https://www.pharmgkb.org/gene/"+l["valueCodeableConcept"]["coding"][0]["code"])

with col2:
    inxData = (pd.DataFrame({
        'Genotype': genotype,
        'Ingredient': medicationName,
        'Ingredient code': medicationCode,
        'Implication': implication,
        'PharmGKB': PharmGKBID,
        'Evidence Level': evidenceLevel}))
    if count == 0:
        st.write('### No PGx data')
    else:
        gb = GridOptionsBuilder.from_dataframe(inxData)
        gb.configure_column("PharmGKB",
                            cellRenderer=JsCode("""
                function (params)
                {return "<a target='_blank' href=" + params.value + ">" + params.value + "</a>"}
                """).js_code
                            )
        go = gb.build()
        AgGrid(inxData, gridOptions=go, allow_unsafe_jscode=True, editable=True, enable_enterprise_modules=True, update_mode="value_changed")
    inxData_t = inxData.T
    st.download_button("Download table (json)", inxData_t.to_json())
    st.download_button("Download table (csv)", inxData.to_csv())

# 🧬
with col1:
    medList = []
    inx = ""
    for i in getMedicationList(subject):
        for j in medicationCode:
            if str(j) == str(i["ing"]):
                inx = '<a target="_blank" href="https://api.pharmgkb.org/v1/infobutton?mainSearchCriteria.v.c='+str(i["ing"])+'">'+'🧬'+'</a>'
                break
            else:
                inx = ""
        medList.append({"-": inx, "Med List": i["product"]})
    data = (pd.DataFrame(medList))
    gb = GridOptionsBuilder.from_dataframe(data)
    gb.configure_column("-",
                        cellRenderer=JsCode("""
             function (params)
             {return params.value}
             """).js_code,
                        )
    go = gb.build()
    AgGrid(data, gridOptions=go, allow_unsafe_jscode=True, enable_enterprise_modules=True, update_mode="value_changed")
