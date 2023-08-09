import { MouseEventHandler, useCallback, useState } from "react";

const dataHeaders = [{
  mnvSPDI: "",
  molecImpact: "",
  snvSPDIs: [""]
}];

type Data = typeof dataHeaders;

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

  const sortedData = tableData.sort((a, b) => {
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

function MNVTable({ data }: { data: Data }) {
  console.log("In sortable table")
  console.log(data)

  const [sortKey, setSortKey] = useState<SortKeys>("mnvSPDI");
  const [sortOrder, setSortOrder] = useState<SortOrder>("ascn");

  const headers: { key: SortKeys; label: string }[] = [
    { key: "mnvSPDI", label: "MNV SPDI" },
    { key: "molecImpact", label: "Molecular Impact" },
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
            <tr key={variant.mnvSPDI}>
              <td>{variant.molecImpact}</td>
            </tr>
          );
        })}
      </tbody>
    </table>
  );
}

export default MNVTable;
