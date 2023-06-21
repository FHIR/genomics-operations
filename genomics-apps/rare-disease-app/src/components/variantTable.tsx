import { useCallback, useState } from "react";

const data = [{ 'id': 1234, 'molecular impact': 'HIGH', 'pathogenicity': 'benign' },
{ 'id': 2345, 'molecular impact': '', 'pathogenicity': 'likely pathogenic' },
{ 'id': 3456, 'molecular impact': 'LOW' },
{ 'id': 4567, 'molecular impact': 'MED', 'pathogenicity': 'benign' },];

const molecImpactDict = { 'MODIFIER': 0, 'LOW': 1, 'MODERATE': 2, 'HIGH': 3 }

type Data = typeof data;

type SortKeys = keyof Data[0];

type SortOrder = 'asc' | 'desc';

function molecImpactSort(a: Data, b: Data) {
    var molecImpactA
    if ('molecular impact' in a) {
        molecImpactA = a['molecular impact']
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

function sortVariables({ tableData, sortKey, reverse }: {
    tableData: Data,
    sortKey: SortKeys,
    reverse: boolean
}) {

    // Not necessary, since automatically sorts by variation id by default
    if (!sortKey) {
        return tableData
    }

    const sortedData = data.sort((a, b) => {
        if (sortKey === 'molecular impact') {
            return molecImpactSort(a, b)
        }

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

function variantTable({ data }: { data: Data }) {

    const [sortKey, setSortKey] = useState<SortKeys>('id');
    const [sortOrder, setSortOrder] = useState<SortOrder>('asc');

    const sortedVariables = useCallback(() =>
        sortVariables({ tableData: data, sortKey, reverse: sortOrder === 'desc' }),
        [data, sortKey, sortOrder])

}
