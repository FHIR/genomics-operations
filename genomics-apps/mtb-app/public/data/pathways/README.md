# Pathway diagrams

The 10 canonical oncogenic signaling pathways from the TCGA PanCancer Atlas, stored as data for the MTB app's Pathways view.

**Source:** Sanchez-Vega F, et al. Oncogenic Signaling Pathways in The Cancer Genome Atlas. *Cell*. 2018;173(2):321-337.e10. doi:10.1016/j.cell.2018.03.035

- Layout, groupings and interactions are redrawn from **Figure 2** (curated pathway templates), as drawn in the figure.
- Gene membership and oncogene / tumor suppressor (OG / TSG) labels come from **Table S3**.

`index.json` lists the pathways. Each pathway file contains:

| Field | Meaning |
|---|---|
| `id`, `name`, `version`, `source` | Identity and provenance |
| `size` | Drawing area (`width`, `height`) in diagram units |
| `nodes` | Boxes. `type` is `gene` or `process` (outputs, stimuli and states such as "Cell growth"). `x`, `y` are the box center; `width`, `height` its size. Gene nodes list their HGNC symbols in `genes`; a combined box such as "PIK3CA/B" lists each gene. |
| `groups` | Boxes around related genes (complexes, families). `x`, `y` are the top-left corner. `members` are node ids; a node can belong to two groups (mTORC1 / mTORC2 share MTOR). `labelPosition` is `top` (default), `bottom-left` or `bottom-right`. |
| `edges` | `from` / `to` are node or group ids. `type` is `activates` (arrow), `inhibits` (flat bar) or `binds` (plain line: binding / part of complex). |
| `genes` | Every Table S3 gene for the pathway: `symbol`, `role` (`OG`, `TSG` or null), and `inDiagram`. Only genes with `inDiagram: true` are drawn and are used for ranking pathways against a patient's variants. |

In the app, a gene box with a patient variant is colored by its `role`: pink for oncogene, blue for tumor suppressor.
