import React, { MouseEventHandler, useCallback, useState } from "react";

function PatientInfoForm({ callback }: {
    callback: (formData: {
        patientID: string;
        geneList: Array<string>;
        addAnnFlag: boolean;
    }) => void
}) {

    const [patientID, setPatientID] = useState('');
    const [geneList, setGeneList] = useState<string[]>([]);
    const [geneString, setGeneString] = useState('BRCA2,APC');
    const [addAnnFlag, setAddAnnFlag] = useState(false);

    const handlePatientIDChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        setPatientID(event.target.value);
    };

    const handleGeneListChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        setGeneString(event.target.value);
        let geneList = event.target.value.split(",");
        setGeneList(geneList);
    };

    const handleAddAnnFlagChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        setAddAnnFlag(event.target.checked);
    };

    const handleSubmit = (event: React.FormEvent<HTMLFormElement>) => {
        event.preventDefault();
        let formData = { patientID: patientID, geneList: geneList, addAnnFlag: addAnnFlag };
        callback(formData);
    }

    return (
        <form onSubmit={handleSubmit}>
            <label htmlFor="patientID">Enter Patient ID:</label>
            <input type="text" id="patientID" name="patientID" value={patientID} onChange={handlePatientIDChange} />

            <label htmlFor="geneString">Enter Gene List, with genes separated by commas:</label>
            <input type="text" id="geneString" name="geneString" value={geneString} onChange={handleGeneListChange} />

            <label htmlFor="addAnnFlag">Compute Additional Annotations:</label>
            <input type="checkbox" name="addAnnFlag" checked={addAnnFlag} onChange={handleAddAnnFlagChange} />

            <button type="submit">Submit</button>
        </form>
    )
}
export default PatientInfoForm
