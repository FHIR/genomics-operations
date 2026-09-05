import { Variant } from '@/types/variants';
import { DxImplication } from '@/services/dxService';
import { ProcessedTxImplication } from '@/services/txService';
import { MolecularConsequence } from '@/services/mcService';
import DxImplicationCell from './DxImplicationCell';
import MolecularConsequenceCell from './MolecularConsequenceCell';
import TxImplicationCell from './TxImplicationCell';
import VariantGroupRow from './VariantGroupRow';
import _ from 'lodash';

interface ResultsTableProps {
  results: Variant[];
}

export default function ResultsTable({ results }: ResultsTableProps) {
  // Function to render Dx Implications content
  const renderDxImplications = (implications?: DxImplication[]) => {
    return <DxImplicationCell implications={implications} />;
  };

  // Function to render Molecular Consequences content
  const renderMolecularConsequences = (consequences?: MolecularConsequence[]) => {
    return <MolecularConsequenceCell consequences={consequences} />;
  };

  // Function to render Tx Implications content
  const renderTxImplications = (implications?: ProcessedTxImplication[]) => {
    return <TxImplicationCell implications={implications} />;
  };

  return (
    <div className="bg-gray-100 p-6 rounded-lg shadow-sm w-full min-w-[1400px] -ml-64">
      <table className="w-full border-collapse table-fixed min-w-[1300px]">
        <colgroup>
          <col style={{width: '200px'}} />
          <col style={{width: '410px'}} />
          <col style={{width: '280px'}} />
          <col style={{width: '235px'}} />
          <col style={{width: '240px'}} />
        </colgroup>
        <thead>
          <tr className="p-2  left bg-gray-200">
            <th className="p-3 text-left text-gray-700">Range</th>
            <th className="p-3 text-left text-gray-700">Variant</th>
            <th className="p-3 text-left text-gray-700">Molecular Consequences</th>
            <th className="p-3 text-left text-gray-700">Dx Implications</th>
            <th className="p-3 text-left text-gray-700">Tx Implications</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(_.groupBy(results, 'range')).map(([range, variantList]) => (
            <VariantGroupRow
              key={range}
              range={range}
              variants={variantList}
              renderDxImplications={renderDxImplications}
              renderTxImplications={renderTxImplications}
              renderMolecularConsequences={renderMolecularConsequences}
            />
          ))}
          {results.length === 0 && (
            <tr>
              <td colSpan={5} className="text-center p-4 text-gray-500">No results found</td>
            </tr>
          )}
        </tbody>
      </table>
    </div>
  );
}