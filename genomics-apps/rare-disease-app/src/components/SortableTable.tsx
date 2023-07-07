import { MouseEventHandler, useCallback, useState } from "react";

const data = [{
  spdi: "",
  dnaChangeType: "",
  sourceClass: "",
  allelicState: "",
  molecImpact: "",
  alleleFreq: NaN,
}];

type Data = typeof data;

type SortKeys = keyof Data[0];

type SortOrder = "ascn" | "desc";

function sortData({
  tableData,
  sortKey,
  reverse,
}: {
  tableData: Data;
  sortKey: SortKeys;
  reverse: boolean;
}) {
  if (!sortKey) return tableData;

  const sortedData = data.sort((a, b) => {
    return a[sortKey] > b[sortKey] ? 1 : -1;
  });

  if (reverse) {
    return sortedData.reverse();
  }

  return sortedData;
}

function SortButton({
  sortOrder,
  columnKey,
  sortKey,
  onClick,
}: {
  sortOrder: SortOrder;
  columnKey: SortKeys;
  sortKey: SortKeys;
  onClick: MouseEventHandler<HTMLButtonElement>;
}) {
  return (
    <button
      onClick={onClick}
      className={`${sortKey === columnKey && sortOrder === "desc"
        ? "sort-button sort-reverse"
        : "sort-button"
        }`}
    >
      ▲
    </button>
  );
}

function SortableTable({ data }: { data: Data }) {
  const [sortKey, setSortKey] = useState<SortKeys>("spdi");
  const [sortOrder, setSortOrder] = useState<SortOrder>("ascn");

  const headers: { key: SortKeys; label: string }[] = [
    { key: "spdi", label: "SPDI" },
    { key: "molecImpact", label: "Molecular Impact" },
    { key: "dnaChangeType", label: "DNA Change Type" },
    { key: "sourceClass", label: "Source Class" },
    { key: "allelicState", label: "Allelic State" },
    { key: "alleleFreq", label: "Allele Frequency" },
  ];

  const sortedData = useCallback(
    () => sortData({ tableData: data, sortKey, reverse: sortOrder === "desc" }),
    [data, sortKey, sortOrder]
  );

  function changeSort(key: SortKeys) {
    setSortOrder(sortOrder === "ascn" ? "desc" : "ascn");

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
        {sortedData().map((variant) => {
          return (
            <tr key={variant.spdi}>
              <td>{variant.spdi}</td>
              <td>{variant.molecImpact}</td>
              <td>{variant.dnaChangeType}</td>
              <td>{variant.alleleFreq}</td>
              <td>{variant.allelicState}</td>
              <td>{variant.sourceClass}</td>
            </tr>
          );
        })}
      </tbody>
    </table>
  );
}

export default SortableTable;
