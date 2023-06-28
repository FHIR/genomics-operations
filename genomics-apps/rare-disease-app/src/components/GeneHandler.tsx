import React, { MouseEventHandler, useCallback, useState } from "react";
import axios from "axios"

const baseURLFSV = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='
const baseURLFC = 'https://fhir-gen-ops.herokuapp.com/utilities/get-feature-coordinates?gene='

const config = {
    headers: {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "GET,PUT,POST,DELETE,PATCH,OPTIONS"
    }
}

function GeneHandler({ patientID, gene, addAnnFlag }:
    {
        patientID: string,
        gene: string,
        addAnnFlag: boolean
    }) {

    console.log("In gene handler")
    console.log("Gene:" + gene)

    var response;

    var range: string;

    console.log("Before API call")
    React.useEffect(() => {
        axios.get(baseURLFC + gene, config).then((response) => {
            console.log(response)
            range = (response.data);
        });
    }, []);
    console.log("After API call")

    React.useEffect(() => {
        if (range) {
            axios.get(baseURLFSV + patientID + '&ranges=' + range + '&includeVariants=true').then((response) => {

            });
        }
    }, []);


    if (patientID === "" || gene === "") return null;


    return null
}
export default GeneHandler
