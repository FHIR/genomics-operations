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
  const patientID = useRef<string>("")
  const [geneList, setGeneList] = useState<string[]>([])
  const addAnnFlag = useRef<boolean>(false)

  var geneTableData: Array<VariantRow> = []

  const handleForm = (formData: { patientID: string, geneList: Array<string>, addAnnFlag: boolean }) => {
    // Update these state variables from within the form component
    patientID.current = formData.patientID
    setGeneList(formData.geneList)
    addAnnFlag.current = formData.addAnnFlag

  }

  const getGeneData = (geneData: Array<VariantRow>) => {
    // Update these state variables from within the form component
    geneTableData = geneData
  }

  return (
    <div className="App">
      <div id="patientInfoForm">
        <PatientInfoForm callback={handleForm} />
      </div>
      {geneList.map(function (gene) {
        console.log(patientID)
        console.log(geneList)
        return <GeneButton patientID={patientID.current} gene={gene} addAnnFlag={addAnnFlag.current} callback={getGeneData} />
      })}

      {/* <SortableTable data={geneTableData} /> */}
    </div>
  );
}

export default App;
