// import React from 'react';
import { MouseEventHandler, useCallback, useState } from "react";
import logo from './logo.svg';
import './App.css';
import SortableTable from "./components/SortableTable";
import PatientInfoForm from "./components/PatientInfoForm";
import GeneHandler from "./components/GeneHandler";
import { ListFormat } from "typescript";
import axios from "axios";

const data = [{ 'id': 1234, 'molecular_impact': 'HIGH', 'pathogenicity': 'benign' },
{ 'id': 2345, 'molecular_impact': '', 'pathogenicity': 'likely pathogenic' },
{ 'id': 3456, 'molecular_impact': 'LOW', 'pathogenicity': 'benign' },
{ 'id': 4567, 'molecular_impact': 'MED', 'pathogenicity': 'benign' },];

function App() {
  const [patientID, setPatientID] = useState("");
  const [geneList, setGeneList] = useState<string[]>([]);
  const [addAnnFlag, setAddAnnFlag] = useState(false);

  const handleForm = (formData: { patientID: string, geneList: Array<string>, addAnnFlag: boolean }) => {
    // Update these state variables from within the form component
    setPatientID(formData.patientID)
    setGeneList(formData.geneList)
    setAddAnnFlag(formData.addAnnFlag)
  }

  return (
    <div className="App">
      <div id="patientInfoForm">
        <PatientInfoForm callback={handleForm} />
        <p>{patientID}</p>
      </div>
      {geneList.map(function (gene) {
        console.log(geneList)
        return <GeneHandler patientID={patientID} gene={gene} addAnnFlag={addAnnFlag} />
      })}

      <SortableTable data={data} />
    </div>
  );
}

export default App;
