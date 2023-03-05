import streamlit as st
import pandas as pd
import requests
import json
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode

st.set_page_config(
    page_title="Simple App",
    page_icon="random",
    layout="wide",
    initial_sidebar_state="collapsed"
)


@st.cache
def findSubjectVariants(subject, range):
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='+subject+'&ranges='+range+'&includeVariants=true'
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)


st.title("Simple App")
st.markdown("This is a simple app illustrating the use of [FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html). Enter patient and genomic range in the sidebar and click 'run'. Raw output from find-subject-variants is shown in left column, and tabular output is shown in right column.")

with st.sidebar:
    subject = st.text_input("Enter patient ID", value="HG00403")
    range = st.text_input("Enter range (zero-based RefSeq:Integer-range)", value="NC_000019.9:11200138-11244496")

col1, col2 = st.columns([2, 3])

if st.sidebar.button("Run"):
    response = findSubjectVariants(subject, range)
    col1.write(response.json())
    if response.status_code == 200:
        variantList = []
        for i in response.json()["parameter"]:
            for j in i["part"]:
                if j["name"] == "variant":
                    sourceClass, allelicState, SPDI, alleleFreq = "", "", "", ""
                    for k in j["resource"]["component"]:
                        if k["code"]["coding"][0]["code"] == "48002-0":
                            sourceClass = k["valueCodeableConcept"]["coding"][0]["display"]
                        elif k["code"]["coding"][0]["code"] == "53034-5":
                            allelicState = k["valueCodeableConcept"]["coding"][0]["display"]
                        elif k["code"]["coding"][0]["code"] == "81252-9":
                            SPDI = k["valueCodeableConcept"]["coding"][0]["display"]
                        elif k["code"]["coding"][0]["code"] == "81258-6":
                            alleleFreq = k["valueQuantity"]["value"]
                    variantList.append({"SPDI": SPDI, "Source Class": sourceClass, "Allelic State": allelicState, "Allele Frequeny": alleleFreq})
        data = (pd.DataFrame(variantList))
        with col2:
            AgGrid(data, enable_enterprise_modules=True, update_mode="value_changed", allow_unsafe_jscode=True)
            st.download_button("Download table (json)", data.T.to_json(), mime="application/json")
            st.download_button("Download table (csv)", data.to_csv(), mime="text/csv")
