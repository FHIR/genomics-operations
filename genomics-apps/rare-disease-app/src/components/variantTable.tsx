import { MouseEventHandler, useCallback, useState } from "react";

const data = [{ 'id': 1234, 'molecular_impact': 'HIGH', 'pathogenicity': 'benign' },
{ 'id': 2345, 'molecular_impact': '', 'pathogenicity': 'likely pathogenic' },
{ 'id': 3456, 'molecular_impact': 'LOW', 'pathogenicity': 'benign' },
{ 'id': 4567, 'molecular_impact': 'MED', 'pathogenicity': 'benign' },];

const molecImpactDict = { 'MODIFIER': 0, 'LOW': 1, 'MODERATE': 2, 'HIGH': 3 }

type Data = typeof data;

type SortKeys = keyof Data[0];

type SortOrder = 'asc' | 'desc';

const headers: { key: SortKeys; label: string }[] = [
    { key: 'id', label: 'ID' },
    { key: 'molecular_impact', label: 'Molecular Impact' },
    { key: 'pathogenicity', label: 'Pathogenicity' },];

// Use if need to list molecular impacts as LOW-HIGH instead of numbers
/*
function molecImpactSort(a: Data, b: Data) {
    var molecImpactA
    if ('molecular impact' in a) {
        molecImpactA = a['molecular impact' as keyof typeof a]
    } else {
        return -1
    }

    var molecImpactB
    if ('molecular impact' in b) {
        molecImpactB = b['molecular impact']
    } else {
        return 1
    }

    return molecImpactDict[molecImpactA] > molecImpactDict[molecImpactB]

}
*/

function sortVariables({ tableData, sortKey, reverse }: {
    tableData: Data,
    sortKey: SortKeys,
    reverse: boolean
}) {

    // Not necessary, since automatically sorts by variation id by default
    if (!sortKey) {
        return tableData
    }

    const sortedData = tableData.sort((a, b) => {
        if (!(sortKey in a) || a[sortKey] === undefined) {
            return -1
        } else if (!(sortKey in b) || b[sortKey] === undefined) {
            return 1
        }
        return (a[sortKey] || '') > (b[sortKey] || '') ? 1 : -1
    }
    )

    if (reverse) {
        return sortedData.reverse
    }

    return sortedData
}

function SortButton({ sortOrder, columnKey, sortKey, onClick }:
    {
        sortOrder: SortOrder
        columnKey: SortKeys
        sortKey: SortKeys
        onClick: MouseEventHandler<HTMLButtonElement>
    }) {
    return <button onClick={onClick}>^</button>
}

function variantTable({ data }: { data: Data }) {

    const [sortKey, setSortKey] = useState<SortKeys>('id');
    const [sortOrder, setSortOrder] = useState<SortOrder>('asc');

    const sortedVariables = useCallback(() =>
        sortVariables({ tableData: data, sortKey, reverse: sortOrder === 'desc' }),
        [data, sortKey, sortOrder])

    function changeSort(key: SortKeys) {
        setSortOrder(sortOrder === 'asc' ? 'desc' : 'asc');

        setSortKey(key);
    }

    return (
        <table>
            <thead>
                <tr>
                    {headers.map((row) => {
                        return (
                            <td key={row.key}>
                                {row.label}{" "}
                                <SortButton
                                    columnKey={row.key}
                                    onClick={() => changeSort(row.key)}
                                    {...{
                                        sortOrder,
                                        sortKey,
                                    }}
                                />
                            </td>
                        );
                    })}
                </tr>
            </thead>

            <tbody>
                {sortedVariables().map((variant) => {
                    return (
                        <tr key={variant.id}>
                            <td>{variant.id}</td>
                            <td>{variant.molecular_impact}</td>
                            <td>{variant.pathogenicity}</td>
                        </tr>
                    );
                })}
            </tbody>
        </table>
    )

}

export default variantTable
