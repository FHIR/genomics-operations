// import React from 'react';
import { useEffect, useState, useRef } from "react";
import logo from './logo.svg';
import './App.css';
import SortableTable from "./components/SortableTable";
import PatientInfoForm from "./components/PatientInfoForm";
import GeneButton from "./components/GeneButton";
import { ListFormat } from "typescript";
import axios from "axios";

const data = [{ 'id': 1234, 'molecular_impact': 'HIGH', 'pathogenicity': 'benign' },
{ 'id': 2345, 'molecular_impact': '', 'pathogenicity': 'likely pathogenic' },
{ 'id': 3456, 'molecular_impact': 'LOW', 'pathogenicity': 'benign' },
{ 'id': 4567, 'molecular_impact': 'MED', 'pathogenicity': 'benign' },];

type VariantRow = {
  spdi: string,
  dnaChangeType: string,
  sourceClass: string,
  allelicState: string,
  molecImpact: string,
  alleleFreq: number,
}

function App() {
  const [geneButtons, setGeneButtons] = useState<{ [key: string]: (Array<VariantRow>) }>({})

  var geneListData: Array<{ geneName: string, geneData?: Array<VariantRow> }> = []

  // const handleForm = (formData: { patientID: string, geneList: Array<string>, addAnnFlag: boolean }) => {
  //   // Update these state variables from within the form component
  //   patientID.current = formData.patientID
  //   setGeneList(formData.geneList)
  //   addAnnFlag.current = formData.addAnnFlag

  // }

  const getGeneData = (geneData: { geneName: string, geneData: Array<VariantRow> }) => {
    // Update state variable from within the form component
    // if (geneData.geneName in geneButtons) {
    setGeneButtons((prevGeneButtons) => {
      return ({
        ...prevGeneButtons,
        geneName: geneData.geneData
      })
    })
    // }
  }

  return (
    <div className="App">
      <div id="patientInfoForm">
        <PatientInfoForm callback={getGeneData} />
      </div>
      {() => {
        if (geneButtons !== undefined && geneButtons !==) {
          geneButtons.keys().map(function (gene) {
            console.log(geneList)
          })
        }
      }}

      <SortableTable data={geneTableData} />
    </div>
  );
}

export default App;
