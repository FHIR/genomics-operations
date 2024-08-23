import pandas as pd
import requests
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder

# Gene ranges and related disease using build37 coordinates
gene_ranges = {
    "APC": {"range": "NC_000005.9:112043194-112181936",
            "disease": "Adenomatous polyposis coli"},
    "MYH11": {"range": "NC_000016.9:15796991-15950885",
              "disease": "Adenomatous polyposis coli"},
    "ACTA2": {"range": "NC_000010.10:90694830-90751154",
              "disease": "Aortic aneurysm, familial thoracic 6"},
    "TMEM43": {"range": "NC_000003.11:14166551-14185180",
               "disease": "Arrhythmogenic right ventricular cardiomyopathy, type 5"},
    "DSP": {"range": "NC_000006.11:7541903-7586947",
            "disease": "Arrhythmogenic right ventricular cardiomyopathy,type 8"},
    "PKP2": {"range": "NC_000012.11:32943688-33049711",
             "disease": "Arrhythmogenic right ventricular cardiomyopathy,\
                         type 9"},
    "DSG2": {"range": "NC_000018.9:29078139-29128971",
             "disease": "Arrhythmogenic right ventricular cardiomyopathy,\
                        type 10"},
    "DSC2": {"range": "NC_000018.9:28638805-28682384",
             "disease": "Arrhythmogenic right ventricular cardiomyopathy,\
                         type 11"},
    "BTD": {"range": "NC_000003.11:15642858-15764023",
            "disease": "Biotinidase deficiency"},
    "BRCA1": {"range": "NC_000017.10:41196311-41277381",
              "disease": "Breast-ovarian cancer, familial 1"},
    "BRCA2": {"range": "NC_000013.10:32889644-32974405",
              "disease": "Breast-ovarian cancer, familial 2"},
    "SCN5A": {"range": "NC_000003.11:38589552-38691178",
              "disease": "Brugada syndrome 1"},
    "RYR2": {"range": "NC_000001.10:237205483-237997288",
             "disease": "Catecholaminergic polymorphic ventricular \
                         tachycardia 1"},
    "CASQ2": {"range": "NC_000001.10:116242641-116311335",
              "disease": "Catecholaminergic polymorphic ventricular \
                          tachycardia 2"},
    "CALM1": {"range": "NC_000014.8:90862845-90874612",
              "disease": "Catecholaminergic polymorphic ventricular \
                          tachycardia 4"},
    "TRDN": {"range": "NC_000006.11:123537483-123958095",
             "disease": "Catecholaminergic polymorphic ventricular \
                         tachycardia 5"},
    "FLNC": {"range": "NC_000007.13:128470459-128499326",
             "disease": "Dilated cardiomyopathy"},
    "LMNA": {"range": "NC_000001.10:156052363-156109872",
             "disease": "Dilated cardiomyopathy 1A"},
    "TNNT2": {"range": "NC_000001.10:201328141-201346808",
              "disease": "Dilated cardiomyopathy 1D"},
    "DES": {"range": "NC_000002.11:220283098-220291456",
            "disease": "Dilated cardiomyopathy 1I"},
    "MYH7": {"range": "NC_000014.8:23881948-23904869",
             "disease": "Dilated cardiomyopathy 1S"},
    "TNNC1": {"range": "NC_000003.11:52485115-52488057",
              "disease": "Dilated cardiomyopathy 1Z"},
    "RBM20": {"range": "NC_000010.10:112404093-112599226",
              "disease": "Dilated cardiomyopathy 1DD"},
    "BAG3": {"range": "NC_000010.10:121410891-121437331",
             "disease": "Dilated cardiomyopathy 1HH"},
    "TTN": {"range": "NC_000002.11:179390715-179672150",
            "disease": "Dilated cardiomyopathy (truncating variants only)"},
    "COL3A1": {"range": "NC_000002.11:189839098-189877472",
               "disease": "Ehlers-Danlos syndrome, type 4"},
    "GLA": {"range": "NC_000023.10:100652790-100662913",
            "disease": "Fabry's disease"},
    "LDLR": {"range": "NC_000019.9:11200138-11244496",
             "disease": "Familial hypercholesterolemia 1"},
    "APOB": {"range": "NC_000002.11:21224300-21266945",
             "disease": "Familial hypercholesterolemia 2"},
    "TPM1": {"range": "NC_000015.9:63334945-63364114",
             "disease": "Familial hypertrophic cardiomyopathy 3"},
    "MYBPC3": {"range": "NC_000011.9:47352956-47374253",
               "disease": "Familial hypertrophic cardiomyopathy 4"},
    "PRKAG2": {"range": "NC_000007.13:151253212-151574200",
               "disease": "Familial hypertrophic cardiomyopathy 6"},
    "TNNI3": {"range": "NC_000019.9:55663134-55669100",
              "disease": "Familial hypertrophic cardiomyopathy 7"},
    "MYL3": {"range": "NC_000003.11:46899361-46904934",
             "disease": "Familial hypertrophic cardiomyopathy 8"},
    "MYL2": {"range": "NC_000012.11:111348648-111358383",
             "disease": "Familial hypertrophic cardiomyopathy 10"},
    "ACTC1": {"range": "NC_000015.9:35082430-35087750",
              "disease": "Familial hypertrophic cardiomyopathy 11"},
    "RET": {"range": "NC_000010.10:43572516-43625799",
            "disease": "Familial medullary thyroid carcinoma"},
    "PALB2": {"range": "NC_000016.9:23614485-23652631",
              "disease": "Hereditary breast cancer"},
    "HFE": {"range": "NC_000006.11:26087656-26098571", "disease":
            "Hereditary hemochromatosis (c.845G>A; p.C282Y homozygotes only)"},
    "ENG": {"range": "NC_000009.11:130577290-130617052",
            "disease": "Hereditary hemorrhagic telangiectasia type 1"},
    "ACVRL1": {"range": "NC_000012.11:52301287-52317145",
               "disease": "Hereditary hemorrhagic telangiectasia type 2"},
    "SDHD": {"range": "NC_000011.9:111957596-111966518",
             "disease": "Hereditary paraganglioma-pheochromocytoma syndrome"},
    "SDHB": {"range": "NC_000001.10:17345216-17380527",
             "disease": "Hereditary paraganglioma-pheochromocytoma syndrome"},
    "TTR": {"range": "NC_000018.9:29171839-29178784",
            "disease": "Hereditary transthyretin-related amyloidosis"},
    "PCSK9": {"range": "NC_000001.10:55505220-55530525",
              "disease": "Hypercholesterolemia, autosomal dominant, 3"},
    "BMPR1A": {"range": "NC_000010.10:88516375-88687726",
               "disease": "Juvenile polyposis syndrome"},
    "SMAD4": {"range": "NC_000018.9:48556582-48611412",
              "disease": "Juvenile polyposis syndrome "},
    "TP53": {"range": "NC_000017.10:7571738-7590808",
             "disease": "Li-Fraumeni syndrome 1"},
    "TGFBR1": {"range": "NC_000009.11:101867394-101916474",
               "disease": "Loeys-Dietz syndrome type 1A"},
    "TGFBR2": {"range": "NC_000003.11:30647993-30735634",
               "disease": "Loeys-Dietz syndrome type 1B"},
    "SMAD3": {"range": "NC_000015.9:67357939-67487507",
              "disease": "Loeys-Dietz syndrome type 3"},
    "KCNQ1": {"range": "NC_000011.9:2466237-2870340",
              "disease": "Long QT syndrome 1"},
    "KCNH2": {"range": "NC_000007.13:150642048-150675409",
              "disease": "Long QT syndrome 2"},
    "CALM2": {"range": "NC_000002.11:47387220-47404075", "disease":
              "Long QT syndrome 15; also associated with catecholaminergic\
                polymorphic ventricular tachycardia"},
    "CALM3": {"range": "NC_000019.9:47104330-47114039", "disease":
              "Long QT syndrome 16; also associated with catecholaminergic\
                polymorphic ventricular tachycardia"},
    "MSH2": {"range": "NC_000002.11:47630205-47710367",
             "disease": "Lynch syndrome 1"},
    "MLH1": {"range": "NC_000003.11:37035008-37092337",
             "disease": "Lynch syndrome 2"},
    "PMS2": {"range": "NC_000007.13:6010555-6048737",
             "disease": "Lynch syndrome 4"},
    "MSH6": {"range": "NC_000002.11:48010283-48034092",
             "disease": "Lynch syndrome 5"},
    "RYR1": {"range": "NC_000019.9:38924330-39078204",
             "disease": "Malignant hyperthermia"},
    "CACNA1S": {"range": "NC_000001.10:201008639-201081554",
                "disease": "Malignant hyperthermia"},
    "FBN1": {"range": "NC_000015.9:48700509-48937906",
             "disease": "Marfan's syndrome"},
    "HNF1A": {"range": "NC_000012.11:121416345-121440315",
              "disease": "Maturity-Onset of Diabetes of the Young"},
    "MEN1": {"range": "NC_000011.9:64570985-64578766",
             "disease": "Multiple endocrine neoplasia, type 1"},
    "MUTYH": {"range": "NC_000001.10:45794913-45806112",
              "disease": "MYH-associated polyposis"},
    "NF2": {"range": "NC_000022.10:29999544-30094589",
            "disease": "Neurofibromatosis, type 2"},
    "OTC": {"range": "NC_000023.10:38211856-38280699",
            "disease": "Ornithine carbamoyltransferase deficiency"},
    "SDHAF2": {"range": "NC_000011.9:61197595-61214205",
               "disease": "Paragangliomas 2"},
    "SDHC": {"range": "NC_000001.10:161284170-161345130",
             "disease": "Paragangliomas 3"},
    "STK11": {"range": "NC_000019.9:1205776-1228430",
              "disease": "Peutz-Jeghers syndrome"},
    "MAX": {"range": "NC_000014.8:65472818-65569413",
            "disease": "Pheochromocytoma"},
    "TMEM127": {"range": "NC_000002.11:96914251-96931735",
                "disease": "Pheochromocytoma"},
    "GAA": {"range": "NC_000017.10:78075379-78093680",
            "disease": "Pompe disease"},
    "PTEN": {"range": "NC_000010.10:89623381-89731687",
             "disease": "PTEN hamartoma tumor syndrome"},
    "RB1": {"range": "NC_000013.10:48877886-49056026",
            "disease": "Retinoblastoma"},
    "RPE65": {"range": "NC_000001.10:68894504-68915637",
              "disease": "RPE65-related retinopathy"},
    "TSC1": {"range": "NC_000009.11:135766735-135820003",
             "disease": "Tuberous sclerosis 1"},
    "TSC2": {"range": "NC_000016.9:2097985-2139492",
             "disease": "Tuberous sclerosis 2"},
    "VHL": {"range": "NC_000003.11:10183461-10195351",
            "disease": "Von Hippel-Lindau syndrome"},
    "WT1": {"range": "NC_000011.9:32409320-32457085",
            "disease": "Wilms' tumor"},
    "ATP7B": {"range": "NC_000013.10:52506804-52585586",
              "disease": "Wilson disease"},
}


# Retrieve patient information
def fetch_patient_info(subject):
    url = f"https://api.logicahealth.org/MTB/open/Patient?identifier={subject}"
    response = requests.get(url)
    data = response.json()
    patient_info = []

    for i in data["entry"]:
        j = i["resource"]
        if j["resourceType"] == "Patient":
            patient_id = j["id"]
            name, dob, gender = "", "", ""
            marital_status, contact, address = "", "", ""
            if j["identifier"][0]["type"]["coding"][0]["display"] == "Medical record number":
                subject_id = j["identifier"][0]["value"]
            name = j["name"][0]["given"][0]
            dob = j["birthDate"]
            gender = j["gender"]
            marital_status = j["maritalStatus"]["text"]
            contact = j["telecom"][0]["value"]
            address = f"{j['address'][0]['line'][0]},\
                {j['address'][0]['city']},\
                {j['address'][0]['state']} - {j['address'][0]['postalCode']},\
                {j['address'][0]['country']}"

            patient_info.append({
                "Patient ID": patient_id,
                "Subject ID": subject_id,
                "Name": name,
                "Gender": gender,
                "Date of Birth": dob,
                "Marital Status": marital_status,
                "Contact": contact,
                "Address": address
            })

    return pd.DataFrame(patient_info)


# Retrieve patient conditions
def fetch_condition(patient_id):
    url = f"https://api.logicahealth.org/MTB/open/Condition?patient={patient_id}"
    response = requests.get(url)
    data = response.json()
    condition_list = []

    if "entry" in data:
        for i in data["entry"]:
            j = i["resource"]
            if j["resourceType"] == "Condition" and \
               j["clinicalStatus"]["coding"][0]["code"] == "active":
                condition = j["code"]["coding"][0]["display"]
                snomed_code = j["code"]["coding"][0]["code"]
                start_date = j["onsetDateTime"].split("T")[0]
                clinical_status = j["clinicalStatus"]["coding"][0]["code"]

                condition_list.append({
                    "Condition": condition,
                    "Snomed Code": snomed_code,
                    "Clinical Status": clinical_status,
                    "Start Date": start_date,
                })

    return pd.DataFrame(condition_list)


# Retrieve medications
def fetch_medication(patient_id):
    url = f"https://api.logicahealth.org/MTB/open/MedicationRequest?patient={patient_id}"
    response = requests.get(url)
    data = response.json()
    medication_list = []

    if "entry" in data:
        for i in data["entry"]:
            j = i["resource"]
            if j["resourceType"] == "MedicationRequest" and j["status"] == "active":
                medication = j["medicationCodeableConcept"]["coding"][0]["display"]
                start_date = j["authoredOn"].split("T")[0]

                medication_list.append({
                    "Medication": medication,
                    "Status": "Active",
                    "Start Date": start_date
                })

    return pd.DataFrame(medication_list)


# Retrieve patient allergies
def fetch_allergy(patient_id):
    url = f"https://api.logicahealth.org/MTB/open/AllergyIntolerance?patient={patient_id}"
    response = requests.get(url)
    data = response.json()
    allergy_list = []

    if "entry" in data:
        for i in data["entry"]:
            j = i["resource"]
            if j["resourceType"] == "AllergyIntolerance":
                substance = j["code"]["coding"][0]["display"]
                reaction = j["reaction"][0]["manifestation"][0]["coding"][0]["display"]

                allergy_list.append({
                    "Substance": substance,
                    "Reaction": reaction
                })

    return pd.DataFrame(allergy_list)


# Retrieve variant information
def fetch_variants(subject, gene):
    url = ("https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?"
           f"subject={subject}&ranges={gene_ranges[gene]['range']}&includeVariants=true&includePhasing=true")
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        variants = []

        for i in data["parameter"]:
            if i["name"] == "variants":
                for j in i["part"]:
                    if j["name"] == "variant":
                        variant_id = j["resource"]["id"]
                        variant_name = None
                        allelic_state = "-"
                        pop_allele_freq = None
                        for comp in j["resource"]["component"]:
                            if comp["code"]["coding"][0]["display"] ==\
                               "Discrete genetic variant":
                                variant_name = comp["valueCodeableConcept"]["coding"][0]["display"]
                            elif comp["code"]["coding"][0]["display"] == "Allelic state":
                                allelic_state = comp["valueCodeableConcept"]["coding"][0]["display"]
                            elif comp["code"]["coding"][0]["display"] ==\
                                    "Population allele frequency":
                                pop_allele_freq = comp["valueQuantity"]["value"]
                        variants.append({
                            "Gene": gene,
                            "Disease Name": gene_ranges[gene]["disease"],
                            "Variant Name": variant_name,
                            "Variant ID": variant_id,
                            "Allelic State": allelic_state,
                            "Population Allele Frequency": pop_allele_freq
                        })

        return variants
    else:
        st.error(f"Failed to fetch variants data: {response.status_code}")
        return []


# Retrieve molecular consequences
def fetch_molecular_consequences(subject, gene):
    url = f"https://fhir-gen-ops.herokuapp.com/subject-operations/phenotype-operations/$find-subject-molecular-consequences?subject={subject}&ranges={gene_ranges[gene]['range']}"

    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()

        if "parameter" not in data:
            return []
        molecular_consq = {}
        priority = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

        for param in data["parameter"]:
            if param["name"] == "consequence":
                j = param["resource"]
                variant_id = j["derivedFrom"][0]["reference"].split("/")[-1]
                impact = j["interpretation"][0]["text"]

                if variant_id not in molecular_consq:
                    molecular_consq[variant_id] = {
                        "Variant ID": variant_id,
                        "Impact": None,
                        "Feature consequence": set(),
                        "Functional effect": None
                    }

                if molecular_consq[variant_id]["Impact"]:
                    current_priority_index = priority.index(molecular_consq[variant_id]["Impact"])
                else:
                    current_priority_index = len(priority)

                # Determine the new priority index
                new_priority_index = priority.index(impact)

                # Update variant information if new impact has higher priority
                if new_priority_index < current_priority_index:
                    molecular_consq[variant_id]["Impact"] = impact
                    molecular_consq[variant_id]["Feature consequence"] = set()
                    molecular_consq[variant_id]["Functional effect"] = None

                # Add feature consequences and functional effect
                # if the new impact has the same or higher priority
                if new_priority_index <= current_priority_index:
                    for component in j["component"]:
                        if component["code"]["coding"][0]["display"] == \
                           "Feature Consequence":
                            feature_consequence = component["valueCodeableConcept"]["coding"][0]["display"]
                            molecular_consq[variant_id]["Feature consequence"].add(feature_consequence)
                        elif component["code"]["coding"][0]["display"] == \
                                "Functional Effect":
                            functional_effect = component["valueCodeableConcept"]["coding"][0]["display"]
                            molecular_consq[variant_id]["Functional effect"] = functional_effect

        final_list = []
        for variant_id, details in molecular_consq.items():
            combined_feature_consq = ";".join(details["Feature consequence"])
            molecular_consq_str = f"{details['Impact']}/{combined_feature_consq}"

            if details["Functional effect"]:
                molecular_consq_str = f"{details['Impact']}/{combined_feature_consq}/{details['Functional effect']}"

            final_list.append({
                "Variant ID": variant_id,
                "Molecular consequence": molecular_consq_str
            })

        return final_list
    else:
        st.error(f"Failed to fetch molecular consequences data: \
                 {response.status_code}")
        return []


# Define function to get level of evidence
def get_level_of_evidence(components):
    evidence_dict = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting classifications": 1,
        "criteria provided, single submitter": 1,
        "no assertion criteria provided": 0,
        "no classification provided": 0,
        "no classification for the individual variant": 0
    }
    level_of_evidence = 0
    for component in components:
        if component["code"]["coding"][0]["display"] == "Level of evidence":
            review_status = component["valueCodeableConcept"]["text"]
            if review_status in evidence_dict:
                level_of_evidence = evidence_dict[review_status]
                break
    return level_of_evidence


# Retrieve pathogenicity based on clinical significance and review status
def fetch_clinical_significance(subject, gene):
    url = f"https://fhir-gen-ops.herokuapp.com/subject-operations/phenotype-operations/$find-subject-dx-implications?subject={subject}&ranges={gene_ranges[gene]['range']}"

    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "parameter" not in data:
            return []
        variants_dict = {}

        for param in data["parameter"]:
            if param["name"] == "implication":
                variant_id = param["resource"]["derivedFrom"][0]["reference"].\
                    split("/")[-1]
                components = param["resource"]["component"]
                clinical_sig = None
                for component in components:
                    if component["code"]["coding"][0]["display"] == \
                       "Genetic variation clinical significance":
                        if "coding" in component["valueCodeableConcept"]:
                            clinical_sig = component["valueCodeableConcept"]["coding"][0]["display"]
                        else:
                            clinical_sig = component["valueCodeableConcept"]["text"]
                        level_of_evidence = get_level_of_evidence(components)
                        clinical_sig_with_evidence = f"{clinical_sig} {level_of_evidence}"
                        if variant_id not in variants_dict:
                            variants_dict[variant_id] = []

                        variants_dict[variant_id].append(clinical_sig_with_evidence)

        # Extract clinical_sig and evidence_list
        clinical_sig_lists = []
        evidence_lists = []
        variant_ids = []

        for variant_id, sig_list in variants_dict.items():
            clinical_sig = []
            evidence = []
            for sig in sig_list:
                parts = sig.rsplit(' ', 1)
                clinical_sig.append(parts[0])
                evidence.append(int(parts[1]))
            clinical_sig_lists.append(clinical_sig)
            evidence_lists.append(evidence)
            variant_ids.append(variant_id)

        # Define groups based on pathogenecity
        concordant_grp1 = ["Benign", "Likely benign"]
        concordant_grp2 = ['Pathogenic', 'Likely pathogenic',
                           'Likely pathogenic, low penetrance',
                           'Pathogenic, low penetrance']

        results = []

        for j in range(0, len(clinical_sig_lists)):
            clinical_sig_list = clinical_sig_lists[j]
            evidence_list = evidence_lists[j]
            clinical_sig_final = []
            variant_id = variant_ids[j]

            max_val = 0
            for i in range(0, len(evidence_list)):
                x = evidence_list[i]
                if x > max_val:
                    max_val = x

            # Collect clinical significance corresponding to max_val
            for i in range(0, len(evidence_list)):
                if evidence_list[i] == max_val:
                    if clinical_sig_list[i] not in clinical_sig_final:
                        clinical_sig_final.append(clinical_sig_list[i])

            clinical_sig_final_str = "|".join(clinical_sig_final)

            if max_val == 4:
                final = f"{clinical_sig_final_str} ; 4 stars"
            elif max_val == 3:
                final = f"{clinical_sig_final_str} ; 3 stars"
            elif max_val == 1:
                if evidence_list.count(max_val) >= 2:
                    if all(clinical_sig in concordant_grp1
                           for clinical_sig in clinical_sig_final):
                        final = f"{clinical_sig_final_str} ; 2 stars"
                    elif all(clinical_sig in concordant_grp2
                             for clinical_sig in clinical_sig_final):
                        final = f"{clinical_sig_final_str} ; 2 stars"
                    else:
                        final = f"{clinical_sig_final_str} ; 1 star"
                else:
                    final = f"{clinical_sig_final_str} ; {max_val} star"
            else:
                final = f"{clinical_sig_final_str} ; {max_val} stars"

            results.append({"Variant ID": variant_id,
                            "Clinical Significance": final})

        return results

    else:
        st.error(f"Failed to fetch clinical significance data: \
                 {response.status_code}")
        return []


# Decorating conditions based on pathogenic variants
def decorate_conditions(condition_df, df_final, selected_genes):
    # Reading the Excel file
    valueset_df = pd.read_excel(r"genomics-apps/data/conditions.xlsx")
    # Convert ValueSetMember to string for comparison
    valueset_df["ValueSetMember"] = valueset_df["ValueSetMember"].astype(str)

    for idx, condition_row in condition_df.iterrows():
        snomed_code = condition_row["Snomed Code"]

        matching_rows = valueset_df[valueset_df["ValueSetMember"] == snomed_code]
        if not matching_rows.empty:
            gene_list = matching_rows["Gene"].str.split(', ')\
                .explode().unique()
            relevant_genes = set(gene_list) & set(selected_genes)
            should_flag = False
            if relevant_genes:
                for gene in relevant_genes:
                    gene_variants = df_final[df_final["Gene"] == gene]
                    for i, variant_row in gene_variants.iterrows():
                        clinical_significance = variant_row.get("Clinical Significance", "")
                        if isinstance(clinical_significance, str):
                            if (("Pathogenic" in clinical_significance or "Likely pathogenic" in clinical_significance)
                                    and ("3 stars" in clinical_significance or "4 stars" in clinical_significance)):
                                should_flag = True
                                break
                    if should_flag:
                        break
            # Add DNA icon to the Condition column only if should_flag is True
            if should_flag:
                condition_df.at[idx, "Condition"] = f"🧬 {condition_df.at[idx, 'Condition']}"
    return condition_df


# Streamlit app
st.set_page_config(
    page_title="Mendelian Screening",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded")

st.title("Detecting Disorders: A Streamlit App")
st.header("Mendelian Disorders Screening", divider="rainbow")
st.markdown('''***An app to find the secondary findings of the disorders.
            This application is a powerful tool designed to assist healthcare
            professionals in screening patients for Mendelian disorders and
            analyzing genetic variants. It provides a comprehensive overview
            of a patient's health status and genetic profile by retrieving and
            displaying essential patient information, active medical
            conditions, current medications, and known allergies. The app's
            core functionality lies in its ability to analyze genetic variants
            for selected genes associated with Mendelian disorders, presenting
            molecular consequences and clinical significance of identified
            variants. By integrating patient health data with genetic analysis,
            this tool aids in the identification and assessment of potential
            genetic disorders, serving as a valuable resource in the realm of
            precision medicine. Conditions associated with pathogenic or
            likely pathogenic variants are flagged with a DNA icon (🧬) for
            easy identification***''')

# Streamlit sidebar for user inputs
st.sidebar.title("Genetic Variant Information")
subject = st.sidebar.text_input("Enter Subject ID")

genes = list(gene_ranges.keys())
selected_genes = st.sidebar.multiselect("Select Genes", genes, default=None)

# Run the app after inputs
if st.sidebar.button("Run"):
    if subject and selected_genes:
        # Row 1: Patient Information
        df = fetch_patient_info(subject)

        st.markdown("""<hr style="height:3px;border:none;color:#333;
        background-color:#333;" /> """, unsafe_allow_html=True)
        st.subheader("Patient Information")

        row1 = st.columns(3)

        if not df.empty:
            with row1[0]:
                st.write("**Patient ID:**", df["Patient ID"][0])
                st.write("**Name:**", df["Name"][0])
                st.write("**Gender:**", df["Gender"][0])

            with row1[1]:
                st.write("**Subject ID:**", df["Subject ID"][0])
                st.write("**Date of Birth:**", df["Date of Birth"][0])
                st.write("**Marital Status:**", df["Marital Status"][0])

            with row1[2]:
                st.write("**Contact:**", df["Contact"][0])
                st.write("**Address:**", df["Address"][0])

        st.markdown("""<hr style="height:3px;border:none;
                    color:#333;background-color:#333;" /> """,
                    unsafe_allow_html=True)

        # Row 2: Condition, Medication, and Allergy
        col1, col2, col3 = st.columns([4, 3, 2])

        with col1:
            condition_df = fetch_condition(df["Patient ID"][0])
            st.subheader("Conditions")
            if condition_df.empty:
                st.write("No conditions found")
            else:
                condition_df_placeholder = st.empty()
                # Display the initial dataframe
                condition_df_placeholder.dataframe(condition_df)

        with col2:
            med_df = fetch_medication(df["Patient ID"][0])
            st.subheader("Medications")
            if med_df.empty:
                st.write("No medications found.")
            else:
                gb = GridOptionsBuilder.from_dataframe(med_df)
                gb.configure_default_column(editable=False, sortable=True,
                                            resizable=True, filterable=True,
                                            filter=True)
                grid_options = gb.build()
                AgGrid(med_df, gridOptions=grid_options,
                       enable_enterprise_modules=True,
                       update_mode="VALUE_CHANGED",
                       allow_unsafe_jscode=True)

        with col3:
            allergy_df = fetch_allergy(df["Patient ID"][0])
            st.subheader("Allergies")
            if allergy_df.empty:
                st.write("No allergies found")
            else:
                gb = GridOptionsBuilder.from_dataframe(allergy_df)
                gb.configure_default_column(editable=False, sortable=True,
                                            resizable=True, filterable=True,
                                            filter=True)
                grid_options = gb.build()
                AgGrid(allergy_df, gridOptions=grid_options,
                       enable_enterprise_modules=True,
                       update_mode="VALUE_CHANGED",
                       allow_unsafe_jscode=True)

        # Row 3: Variants Table
        st.markdown("""<hr style="height:3px;border:none;color:#333;
                    background-color:#333;" /> """, unsafe_allow_html=True)
        st.subheader("Genetic Variant Information")

        st.write(f"Fetching variants for Subject ID: **{subject}** and \
                 Gene(s): **{', '.join(selected_genes)}**...")

        with st.spinner("Fetching data..."):
            all_variants = []
            all_molecular_consequences = []
            all_clinical_significance = []

            for gene in selected_genes:
                variants = fetch_variants(subject, gene)
                molecular_consequence = fetch_molecular_consequences(subject, gene)
                clinical_significance = fetch_clinical_significance(subject, gene)

                if variants:
                    all_variants.extend(variants)
                else:
                    st.warning(f"No variant found for gene {gene}")

                if molecular_consequence:
                    all_molecular_consequences.extend(molecular_consequence)
                else:
                    st.warning(f"No molecular consequence \
                               found for gene {gene}")
                if clinical_significance:
                    all_clinical_significance.extend(clinical_significance)
                else:
                    st.warning(f"No clinical significance \
                               found for gene {gene}")

            if all_variants:
                st.write(f"**Total variants found: {len(all_variants)}**")
                df_variants = pd.DataFrame(all_variants)
                df_molecular_consq = pd.DataFrame(all_molecular_consequences)
                df_clinical_sig = pd.DataFrame(all_clinical_significance)
                # Merge dataframes
                if not df_molecular_consq.empty and \
                   'Variant ID' in df_variants.columns and \
                   'Variant ID' in df_molecular_consq.columns:
                    df_merged = pd.merge(df_variants, df_molecular_consq, how="left", on="Variant ID")
                else:
                    df_merged = df_variants

                if not df_clinical_sig.empty and \
                   'Variant ID' in df_merged.columns and \
                   'Variant ID' in df_clinical_sig.columns:
                    df_final = pd.merge(df_merged, df_clinical_sig, how="left", on="Variant ID")
                else:
                    df_final = df_merged

                # Drop the 'Variant ID' column
                df_final = df_final.drop(columns=["Variant ID"])

                # Configure AgGrid options for the variant dataframe
                gb = GridOptionsBuilder.from_dataframe(df_final)
                gb.configure_pagination()
                gb.configure_default_column(editable=False, sortable=True,
                                            resizable=True, filterable=True,
                                            filter=True)
                grid_options = gb.build()

                st.write("#### Variants with Molecular Consequences \
                         and Pathogenicity")
                grid_response = AgGrid(
                    df_final,
                    gridOptions=grid_options,
                    enable_enterprise_modules=True,
                    update_mode="VALUE_CHANGED",
                    allow_unsafe_jscode=True
                )

                # Update the conditions based on pathogenic variant
                updated_condition_df = decorate_conditions(condition_df, df_final, selected_genes)

                # Drop the 'Snomed Code' column
                updated_condition_df = updated_condition_df.drop(
                    columns=["Snomed Code"])

                # Configure AgGrid options for the updated dataframe
                gb = GridOptionsBuilder.from_dataframe(updated_condition_df)
                gb.configure_default_column(editable=True, sortable=True,
                                            resizable=True, filterable=True,
                                            filter=True)
                grid_options = gb.build()

                # Clear the placeholder
                condition_df_placeholder.empty()

                # Display the updated dataframe with AgGrid in the placeholder
                with condition_df_placeholder.container():
                    AgGrid(updated_condition_df, gridOptions=grid_options,
                           enable_enterprise_modules=True,
                           update_mode="VALUE_CHANGED",
                           allow_unsafe_jscode=True
                           )
            else:
                st.warning("No variants found for the given inputs.")
    else:
        st.warning("Please enter both Subject ID \
                   and select at least one Gene.")
else:
    st.write("Please enter a Subject ID, select Genes, \
             and click 'Run' in the sidebar to start the analysis.")

st.markdown("---")
