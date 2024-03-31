import streamlit as st
import requests
import csv
import statistics

st.set_page_config(
    page_title="Get Polygenic Score",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.title("🧬Get Polygenic Score")
st.markdown("This proof-of-concept app demonstrates the computation of a polygenic score, \
    **using [FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html) to access a person's entire genome**.  \n\
    Select patient and polygenic model in the sidebar, optionally adjust standard deviation threshold, and click 'Run'.")
with st.expander("SIMPLIFICATIONS..."):
    st.markdown("\
        :one: **Wildtype:** We assume wildtype if no variant has been reported at a position within the studied region(s).  \n\
        :two: **Missing Data:** Imputation / statistical approaches for missing data are deferred for now. We currently assume wildtype for polymorphisms outside the studied region(s).  \n\
        :three: **Ancestry:** No ancestry considerations are yet factored in.  \n\
        :four: **Normalization:** Alleles are normalized to [canonical SPDI](https://www.ncbi.nlm.nih.gov/variation/notation/) to enhance search.  \n\
        :five: **Performance:** Proof-of-concept code is not optimized for performance.  \n\
        :six: **Scoring:**  \n\
            :arrow_forward: Loosely follows the methods described by Hao, et al ([PMID: 35437332](https://pubmed.ncbi.nlm.nih.gov/35437332/)).  \n\
            :arrow_forward: *score:* product of germline allele dose (absent=0; heterozygous=1; homozygous=2) x effect_weight.  \n\
            :arrow_forward: *polygenic score (raw):* sum of scores.  \n\
            :arrow_forward: *average risk:* polygenic score (raw) <= selected StDev from the population mean.  \n\
            :arrow_forward: *high risk:* polygenic score (raw) > selected StDev from the population mean.   \n\
            :arrow_forward: *population mean* and *standard deviation* are based on Operations reference implementation population.")

@st.cache
def findSubjectVariants(subject, range):  
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='+subject+'&ranges='+range+'&includeVariants=true'
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)

@st.cache
def findSubjectSpecificVariants(subject, variant):
    url = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-specific-variants?subject='+subject+'&variants='+variant+'&includeVariants=true'
    # url = 'https://fhir-gen-ops-dev-ca42373833b6.herokuapp.com/subject-operations/genotype-operations/$find-subject-specific-variants?subject='+subject+'&variants='+variant+'&includeVariants=true'
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)

@st.cache
def getPolygenicScoreMetadata(pgs_id):
    url = 'https://www.pgscatalog.org/rest/score/'+pgs_id
    headers = {'Accept': 'application/json'}
    return requests.get(url, headers=headers)

def populatePolygenicScoreModels():
    with open('genomics-apps/data/polygenicScore.csv') as polygenicScoreFile:
        polygenicScores = csv.reader(polygenicScoreFile, delimiter=',', quotechar='"')
        pgs_id_list = []
        phenotype_list = []
        pgs_name_list = []
        GRCh37_SPDI_list = []
        effect_weight_list = []
        for item in polygenicScores:
            if item[0] != "pgs_id":
                pgs_id_list.append(item[0])
                phenotype_list.append(item[1])
                GRCh37_SPDI_list.append(item[3])
                effect_weight_list.append(item[10])
    polygenicScoreModels = [pgs_id_list, phenotype_list, GRCh37_SPDI_list, effect_weight_list]
    return polygenicScoreModels

def populatePolygenicPopulationStatistics():
    with open('genomics-apps/data/polygenicPopulationStatistics.csv') as polygenicPopulationStatisticsFile:
        populationStatistics = csv.reader(polygenicPopulationStatisticsFile, delimiter=',', quotechar='"')
        pgs_id_list = []
        mean_list = []
        stdev_list = []
        for item in populationStatistics:
            pgs_id_list.append(item[0])
            mean_list.append(item[1])
            stdev_list.append(item[2])
    polygenicPopulationStatistics = [pgs_id_list, mean_list, stdev_list]
    return polygenicPopulationStatistics

def getPolygenicRawScore (subject, polygenicModelID):
    itemNumber = 0
    phenotype = ""
    results_table = []
    total_score = 0
    for item in polygenicScoreModels[0]:
        if item == polygenicModelID:
            phenotype = polygenicScoreModels[1][itemNumber]
            GRCh37_SPDI = polygenicScoreModels[2][itemNumber]
            effect_weight = float(polygenicScoreModels[3][itemNumber])
            if GRCh37_SPDI.split(":")[2]!=GRCh37_SPDI.split(":")[3]: # effect allele is not wildtype
                variant = findSubjectSpecificVariants(subject,GRCh37_SPDI)
                print(f"GRCh37_SPDI: {GRCh37_SPDI}; variant: {variant}") # returns an error if NCBI SPDI services are down
                if variant.json()["parameter"][0]["part"][1]["valueBoolean"]==False: 
                    dose=0
                    score=0
                    total_score = total_score + score
                    results_table.append({"allele":GRCh37_SPDI,"dose":dose,"effect_weight":effect_weight,"score":score})
                else:
                    for j in variant.json()["parameter"][0]["part"][2]["resource"]["component"]:
                        if j["code"]["coding"][0]["code"] == "53034-5":
                            allelicStateCode = j["valueCodeableConcept"]["coding"][0]["code"]                       
                            allelicStateDisplay = j["valueCodeableConcept"]["coding"][0]["display"]
                            if allelicStateCode == "LA6706-1":
                                dose=1
                                score=dose * effect_weight
                                results_table.append({"allele":GRCh37_SPDI,"dose":dose,"effect_weight":effect_weight,"score":score})
                                total_score = total_score + score
                            elif allelicStateCode == "LA6705-3":
                                dose=2
                                score=dose * effect_weight
                                results_table.append({"allele":GRCh37_SPDI,"dose":dose,"effect_weight":effect_weight,"score":score})
                                total_score = total_score + score
                            else:
                                dose="undefined"
                                score="undefined"
                                results_table.append({"allele":GRCh37_SPDI,"dose":"undefined","effect_weight":effect_weight,"score":"undefined"})
            else: # effect allele is wildtype
                refSeq = GRCh37_SPDI.split(":")[0]
                position = GRCh37_SPDI.split(":")[1]
                variant = findSubjectVariants(subject,refSeq+":"+position+"-"+position)
                if variant.json()["parameter"][0]["part"][1]["valueBoolean"]==False: # no variant found, thus homozygous wildtype
                    dose=2
                    score=dose * effect_weight
                    results_table.append({"allele":GRCh37_SPDI,"dose":dose,"effect_weight":effect_weight,"score":score})
                    total_score = total_score + score
                else:
                    for j in variant.json()["parameter"][0]["part"][2]["resource"]["component"]:
                        if j["code"]["coding"][0]["code"] == "53034-5":
                            allelicStateCode = j["valueCodeableConcept"]["coding"][0]["code"]                       
                            allelicStateDisplay = j["valueCodeableConcept"]["coding"][0]["display"]
                            if allelicStateCode == "LA6706-1":
                                dose=1
                                score=dose * effect_weight
                                results_table.append({"allele":GRCh37_SPDI,"dose":dose,"effect_weight":effect_weight,"score":score})
                                total_score = total_score + score
                            elif allelicStateCode == "LA6705-3":
                                dose=0
                                score=dose * effect_weight
                                results_table.append({"allele":GRCh37_SPDI,"dose":dose,"effect_weight":effect_weight,"score":score})
                                total_score = total_score + score
                            else:
                                dose="undefined"
                                score="undefined"
                                results_table.append({"allele":GRCh37_SPDI,"dose":"undefined","effect_weight":effect_weight,"score":"undefined"})
        itemNumber = itemNumber + 1
    polygenicRawScore = {"subject": subject, "polygenicModelID":polygenicModelID, "phenotype": phenotype, "scores": results_table, "rawScore": total_score}
    return polygenicRawScore

def getSummaryStatistics (polygenicModelID):
    score = []
    for subject in ["HG00403","HG00406","HG02657","NA18498","NA18499","NA18870","NA18871","NA19190","NA19210","NA19238","NA19239","NA19240"]:
        rawScore = getPolygenicRawScore(subject, polygenicModelID)["rawScore"]
        score.append(rawScore)
    mean = statistics.mean(score)
    stdev = statistics.stdev(score)
    polygenicSummaryStatistics = {"polygenicModelID": polygenicModelID, "mean": mean, "stdev": stdev}
    return polygenicSummaryStatistics

polygenicScoreModels = populatePolygenicScoreModels()
polygenicPopulationStatistics = populatePolygenicPopulationStatistics()

# ******* This code computes mean and stdev for ALL models - it takes a while to run ******
# ******* Alternatively, call getSummaryStatistics for just the model of interest *********
# for polygenicModelID in set(polygenicScoreModels[0]):
    # print(getSummaryStatistics(polygenicModelID))
# print(getSummaryStatistics("PHECODE15S"))

with st.sidebar:
    subject = st.selectbox("Select patient", ["HG00403","HG00406","HG02657","NA18498","NA18499","NA18870","NA18871","NA19190","NA19210","NA19238","NA19239","NA19240"])
    unique_pgs_id_list = []
    for item in polygenicScoreModels[0]:
        if item not in unique_pgs_id_list:
            unique_pgs_id_list.append(item)
    polygenicModelID = st.selectbox("Select polygenic model", unique_pgs_id_list)
    riskThreshold = st.slider("Select StDev threshold for high risk",max_value=3.0,value=1.0,step=0.1)
    st.image("genomics-apps/data/normalDistribution.png",width=300)

if st.sidebar.button("Run"):
    polygenicRawScore = getPolygenicRawScore(subject, polygenicModelID)
    rawScore = float(polygenicRawScore["rawScore"])
    polygenicModelID = polygenicRawScore["polygenicModelID"]
    phenotype = polygenicRawScore["phenotype"]

    itemCount = 0
    for item in polygenicPopulationStatistics[0]:
        if item == polygenicModelID:
            populationMean = float(polygenicPopulationStatistics[1][itemCount])
            populationStDev = float(polygenicPopulationStatistics[2][itemCount])
        itemCount = itemCount + 1
    polygenicRisk = (rawScore - populationMean) / (populationStDev + 0.0000001)

    st.subheader(f"{polygenicModelID} ({phenotype})")
    if polygenicRisk > riskThreshold:
        st.image("genomics-apps/data/highRisk.png", width=300)
    else:
        st.image("genomics-apps/data/averageRisk.png", width=300)

    st.markdown(f"\
        **Polygenic score (raw):** {rawScore:.4f}  \n\
        **Population mean:** {populationMean:.4f}  \n\
        **Population standard deviation:** {populationStDev:.4f}  \n\
        **{polygenicRisk:.2f}** standard deviations above the mean.")

    with st.expander("Allele-based scores"):
        st.table(data=polygenicRawScore["scores"])

    with st.expander("Polygenic Score Metadata"):
        st.write(getPolygenicScoreMetadata(polygenicModelID).json())
   
    with st.expander("Polygenic Score Model"):
        polygenicModelFile = "genomics-apps/data/" + polygenicModelID + ".txt"
        polygenicModel = open("genomics-apps/data/"+ polygenicModelID + ".txt","r")
        fileContents=""
        for item in polygenicModel:
            fileContents=fileContents+item
        st.text(fileContents)