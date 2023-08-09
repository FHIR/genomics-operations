import { getEventListeners } from "events";
import React from "react";

type VariantRow = {
    spdi: string,
    dnaChangeType: string,
    sourceClass: string,
    allelicState: string,
    molecImpact: string,
    alleleFreq: number,
}

type MNVRow = {
    mnvSPDI: string,
    molecImpact: string,
    snvSPDIs: Array<string>
}

type FSVResponse = {
    parameter: Array<{
        name: string,
        part: Array<FSVParameter>
    }>,
    resourceType: string
}

type FHIRCoding = [{
    code: string,
    system: string,
    display?: string
}]

type FSVVariantResource = {
    category: Array<FHIRCoding>,
    code: FHIRCoding,
    component: Array<{
        code: { coding: FHIRCoding },
        valueCodeableConcept?: { coding: FHIRCoding },
        valueString?: string,
        interpretation?: { text: string },
        valueQuantity?: { code: string, system: string, value: number },
        valueRange?: { low: { value: number } }
    }>,
    id: string,
    meta: { profile: Array<string> },
    resourceType: string,
    status: string,
    subject: { reference: string },
    valueCodeableConcept: { coding: FHIRCoding },
    text?: string,
}

type PhaseRelationshipResource = {
    valueCodeableConcept: { coding: FHIRCoding },
    derivedFrom: Array<{ reference: string }>
}

type FSVParameter = {
    name: string,
    resource?: FSVVariantResource | PhaseRelationshipResource,
    valueString?: string,
    valueBoolean?: boolean,
}

type variant = {
    id: string,
}

const baseURLFSV = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='
const baseURLFC = 'https://fhir-gen-ops.herokuapp.com/utilities/get-feature-coordinates?gene='

const config = {
    headers: {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "GET,PUT,POST,DELETE,PATCH,OPTIONS"
    }
}

function translateComponents(resource: FSVVariantResource) {
    let variantRow: VariantRow = {
        spdi: "",
        dnaChangeType: "",
        sourceClass: "",
        allelicState: "",
        molecImpact: "",
        alleleFreq: NaN,
    }
    resource.component.map((component) => {
        if (component.code.coding[0].code == "48002-0") {
            if (component.valueCodeableConcept && component.valueCodeableConcept.coding[0].display) {
                variantRow.sourceClass = component.valueCodeableConcept.coding[0].display
            }
        } else if (component.code.coding[0].code == "53034-5") {
            if (component.valueCodeableConcept && component.valueCodeableConcept.coding[0].display) {
                variantRow.allelicState = component.valueCodeableConcept.coding[0].display
            }
        } else if (component.code.coding[0].code == "81252-9") {
            if (component.valueCodeableConcept && component.valueCodeableConcept.coding[0].display) {
                variantRow.spdi = component.valueCodeableConcept.coding[0].display
            }
        } else if (component.code.coding[0].code == "92821-8") {
            if (component.valueQuantity && component.valueQuantity.value) {
                variantRow.alleleFreq = component.valueQuantity.value
            }
        } else if (component.code.coding[0].code == "molecular-consequence") {
            if (component.interpretation) {
                variantRow.molecImpact = component.interpretation.text
            }
        }
    })
    return variantRow
}

function parsePhaseRelationship(param: FSVParameter, addAnnFlag: boolean) {
    if (param.resource?.valueCodeableConcept.coding[0].code != 'Cis') {
        return
    }

    let resource = param.resource as PhaseRelationshipResource
    let derivedFromIDs: Array<string> = []

    resource.derivedFrom.map((idEntry) => {
        derivedFromIDs.push(idEntry.reference.split('/')[1])
    })

    return derivedFromIDs
}

async function computeAnnotations(SPDI: string) {
    type SPDISeqInfo = {
        seq_id: string,
        position: number,
        deleted_sequence: string,
        inserted_sequence: string,
        warnings?: [
          {
            reason: string,
            message: string
          }
        ]
    }
    type AllEquivSPDIResponse = {data?: {
        spdis: [SPDISeqInfo, SPDISeqInfo],
        errors: {
          error: {
            code: number,
            message: string,
            errors: [
              {
                reason: string,
                message: string
              }
            ]
          }
        },
        warnings: [
          {
            reason: string,
            message: string
          }
        ]
      }
    }
    type HGVSResponse = {
        data: {
          hgvs: string
          warnings: [
            {
              reason: string,
              message: string
            }
          ]
        }
    }
    type vepResponse = [
        {
            seq_region_name: string,
            transcript_consequences: [
                {
                    impact: string
                }
            ]
        }

    ]
    var url = `/spdi/${SPDI.replace(":", "%3A")}/all_equivalent_contextual`

    var response = await fetch(url)
    var responseJson = await response.json() as AllEquivSPDIResponse
    if (responseJson instanceof Error) {
        console.log('It is an error!');
    }
    else {
        console.log(responseJson);
    }

    let b38SeqID: string = ""
    let b38Pos: number = -1
    if (responseJson.data) {
        b38SeqID = responseJson.data.spdis[1].seq_id
        b38Pos = responseJson.data.spdis[1].position
    }

    if (b38SeqID == "" || b38Pos == -1) {
        console.log("Error retrieving SPDI info")
        return ""
    }

    var SPDIList = SPDI.split(":")
    let b38SPDI = b38SeqID + ":" + b38Pos + ":" + SPDIList[2] + ":" + SPDIList[3]

    var url = 'https://api.ncbi.nlm.nih.gov/variation/v0/spdi/'+ b38SPDI.replace(":", "%3A") + '/hgvs'

    var response = await fetch(url)
    var HGVSResponseJson = await response.json() as HGVSResponse
    if (HGVSResponseJson instanceof Error) {
        console.log('Error getting HGVS');
    }
    else {
        console.log(HGVSResponseJson);
        var HGVS = HGVSResponseJson.data.hgvs
    }

    var response = await fetch(url)
    var vepResponseJson = await response.json() as vepResponse
    if (vepResponseJson instanceof Error) {
        console.log('Error getting vep response');
    }
    else {
        console.log(vepResponseJson);
    }

    let molecImpact = ""

    vepResponseJson[0].transcript_consequences.map((dict) => {
        if ("impact" in dict && molecImpact == "") {
            molecImpact = "Computed: " + dict.impact
        }
    })


    return molecImpact

    // Call variation services all equivalent
    // Take [1] entry
    // Send to hgvs
}

function findMNVs(response: FSVResponse, cisVariantsIDs: Array<Array<string>>, mnvData: Array<{
    mnvSPDI: string,
    molecImpact: string,
    snvSPDIs: Array<string>
}>) {
    cisVariantsIDs.map(async (IDList) => {
        let SPDIList: Array<string> = []
        let SPDIDictList: Array<{ SPDI: string, chromosome: string, position: string, refAllele: string, altAllele: string }> = []
        response.parameter[0].part.map((param) => {
            // Return if not variant profile
            if (param.name != 'variant' || param.resource == null) {
                return
            }
            // Typecast resource as variant resource
            let resource = param.resource as FSVVariantResource

            // Return if id is not part of an identified mnv
            if (!IDList.includes(resource.id)) {
                return
            }

            resource.component.map((component) => {
                // Return if not spdi component
                if (component.code.coding[0].code != '81252-9' || !component.valueCodeableConcept || !component.valueCodeableConcept.coding[0].display) {
                    return
                }

                // Return if already parsed
                if (SPDIList.includes(component.valueCodeableConcept.coding[0].display)) {
                    return
                }

                let SPDIString = component.valueCodeableConcept.coding[0].display
                SPDIList.push(SPDIString)
                let SPDIStringList = SPDIString.split(':')
                SPDIDictList.push({
                    SPDI: SPDIString,
                    chromosome: SPDIStringList[0],
                    position: SPDIStringList[1],
                    refAllele: SPDIStringList[2],
                    altAllele: SPDIStringList[3],
                })
            })
        })

        // Sort the spdis based on position
        SPDIDictList.sort(function (a, b) {
            if (a.position < b.position) return -1;
            if (a.position > b.position) return 1;
            return 0;
        });

        let refAllele = ""
        let altAllele = ""

        SPDIDictList.map((SPDIDict) => {
            refAllele += SPDIDict.refAllele
            altAllele += SPDIDict.altAllele
        })

        let mnvSPDI = SPDIDictList[0].chromosome + ":" + SPDIDictList[0].position +
            ":" + refAllele + ":" + altAllele

        mnvData.push({
            mnvSPDI: mnvSPDI,
            snvSPDIs: [SPDIDictList[0].SPDI, SPDIDictList[1].SPDI],
            molecImpact: await computeAnnotations(mnvSPDI)
        })
    })

    return mnvData
}

function translateFHIRResponse(response: FSVResponse, addAnnFlag: boolean) {
    let geneTable: Array<VariantRow> = []
    let paramArray = response.parameter[0].part;
    let cisVariantsIDs: Array<Array<string>> = []
    let mnvData: Array<MNVRow> = []
    paramArray.map((param) => {
        //If the part is not a variant resource, return
        if (param.name == 'sequencePhaseRelationship' && addAnnFlag) {
            let derivedFromIDs = parsePhaseRelationship(param, addAnnFlag)
            if (derivedFromIDs) {
                cisVariantsIDs.push(derivedFromIDs)
            }
            return
        } else if (param.name != 'variant' || param.resource == null) {
            return
        }

        let variantRow: VariantRow = translateComponents(param.resource as FSVVariantResource)
        geneTable.push(variantRow)
    })

    if (addAnnFlag) {
        findMNVs(response, cisVariantsIDs, mnvData)
    }

    return {geneTable: geneTable, mnvData: mnvData}
}

function getGeneData({ patientID, gene, addAnnFlag, callback }:
    {
        patientID: string,
        gene: string,
        addAnnFlag: boolean,
        callback: (geneData: { geneName: string, geneData: Array<VariantRow>, mnvData: Array<MNVRow> }) => void
    }) {

    // const [APIStatus, setAPIStatus] = useState("red")
    var APIStatus = "red"
    var geneData: Array<VariantRow>
    var mnvData: Array<MNVRow>

    console.log("In gene handler")
    console.log("Gene:" + gene)

    console.log("PatientID: " + patientID)

    if (!patientID) {
        return null
    }

    async function getFHIRResponse() {
        let range = ""
        type FeatureCoordinates = [{
            MANE: string[],
            build37Coordinates: string,
            build38Coordinates: string,
            geneId: string,
            geneLink: string,
            geneSymbol: string,
            transcripts: string[]
        }]

        let url = baseURLFC + gene
        // urlAppender = urlAppender.replaceAll("/", "@")
        // urlAppender = urlAppender.replaceAll("$", "@!abab@!")
        // urlAppender = urlAppender.replaceAll("?", "!")


        // let url = `http://127.0.0.1:5000/${urlAppender}`

        let response: any

        response = await fetch(url)
        var responseJson = await response.json() as FeatureCoordinates
        if (responseJson instanceof Error) {
            console.log('It is an error!');
        }
        else {
            console.log(responseJson);
            // setAPIStatus("yellow")
            APIStatus = 'yellow'
            // setArticleDict(JSON.parse(responseJson));
        }

        range = responseJson[0].build38Coordinates

        if (range == "") {
            return
        }

        url = baseURLFSV + patientID + '&ranges=' + range + '&includeVariants=true' + '&includePhasing=true'
        // urlAppender = urlAppender.replaceAll("/", "@")
        // urlAppender = urlAppender.replaceAll("$", "@!abab@!")
        // urlAppender = urlAppender.replaceAll("?", "!")

        // let url = `http://127.0.0.1:5000/${urlAppender}`
        console.log(url)

        let fsvResponse = await fetch(url)
        var fsvResponseJson = await fsvResponse.json() as FSVResponse
        if (fsvResponseJson instanceof Error) {
            console.log('It is an error!');
        }
        else {
            console.log(fsvResponseJson.parameter[0]);
            let tableData = await translateFHIRResponse(fsvResponseJson, addAnnFlag);
            geneData = tableData.geneTable
            mnvData = tableData.mnvData
            // setAPIStatus("green")
            APIStatus = 'green'
        }

        callback({ geneName: gene, geneData: geneData, mnvData: mnvData })
    }

    getFHIRResponse()
}

export default function PatientInfoForm({ callback }: {
    callback: (geneData: { geneName: string, geneData: Array<VariantRow>, mnvData: Array<MNVRow> }) => void
}) {

    function handleSubmit(event: React.FormEvent<HTMLFormElement>) {
        interface FormDataElements extends HTMLFormControlsCollection {
            patientID: HTMLInputElement;
            geneList: HTMLInputElement;
            addAnnFlag: HTMLInputElement;
        }

        event.preventDefault();
        const elements = event.currentTarget.elements as FormDataElements;

        const data = {
            patientID: elements.patientID.value,
            geneList: elements.geneList.value.split(","),
            addAnnFlag: elements.addAnnFlag.checked,
        };



        console.log(`Here's your data: ${JSON.stringify(data, undefined, 2)}`);
        data.geneList.map((gene) => {
            getGeneData({ patientID: data.patientID, gene: gene, addAnnFlag: data.addAnnFlag, callback: callback })
        })

    }

    return (
        <form onSubmit={handleSubmit}>
            <label htmlFor="patientID">Patient ID</label>
            <input id="patientID" type="text" />

            <label htmlFor="geneList">Gene List</label>
            <input id="geneList" type="text" />

            <input id="addAnnFlag" type="checkbox" />
            <label htmlFor="addAnnFlag">Compute Additional Annotations</label>

            <button type="submit">Submit</button>
        </form>
    );
}
