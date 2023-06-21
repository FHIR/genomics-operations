import React from 'react';
import logo from './logo.svg';
import './App.css';
import SortableTable from "./components/SortableTable";

const data = [{ 'id': 1234, 'molecular_impact': 'HIGH', 'pathogenicity': 'benign' },
{ 'id': 2345, 'molecular_impact': '', 'pathogenicity': 'likely pathogenic' },
{ 'id': 3456, 'molecular_impact': 'LOW', 'pathogenicity': 'benign' },
{ 'id': 4567, 'molecular_impact': 'MED', 'pathogenicity': 'benign' },];

function App() {
  return (
    <div className="App">
      <SortableTable data={data} />
    </div>
  );
}

export default App;
