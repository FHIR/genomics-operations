# Design: `group` parameter for population operations

Status: **design finalized; groups loaded; awaiting go-ahead to code** (no code written yet)
Branch: `add-group-parameter`
Source spec: HL7 Genomics Reporting IG, FHIR-59649 "Specify Denominator for Population Queries"
- Overview: https://build.fhir.org/ig/HL7/genomics-reporting/branches/operations/operations.html
- Example OperationDefinition: https://build.fhir.org/ig/HL7/genomics-reporting/branches/operations/OperationDefinition-find-population-specific-variants.html

## 1. What the IG requires

All seven population operations get a new optional input parameter:

| Name | Use | Card. | Type | Meaning |
| :--- | :--- | :--- | :--- | :--- |
| `group` | in | 0..1 | string | Identifier of a FHIR Group already known to the server. |

When `group` is supplied:
- **numerator** counts only group members that match the query
- **denominator** is the number of members in the group
- **subject** list (when `includePatientList=true`) contains only group members
- an unknown group identifier returns **404 `ERROR: Group not found`**

When `group` is omitted, behavior is unchanged: the cohort is every patient in the repository.

Requesters can't define groups or pass patient lists through the operation. They can only name a group that the server already holds.

Affected operations (all in [app/endpoints.py](app/endpoints.py)):

| Operation | Handler | Result blocks with numerator/denominator |
| :--- | :--- | :--- |
| find-population-specific-variants | `find_population_specific_variants` | 2 (multi-param AND; per-item) |
| find-population-structural-intersecting-variants | `find_population_structural_intersecting_variants` | 1 (per range) |
| find-population-structural-subsuming-variants | `find_population_structural_subsuming_variants` | 1 (per range) |
| find-population-specific-haplotypes | `find_population_specific_haplotypes` | 2 (multi-param AND; per-item) |
| find-population-tx-implications | `find_population_tx_implications` | 4 (variants; haplotypes; treatments; conditions) |
| find-population-dx-implications | `find_population_dx_implications` | 3 (variants; haplotypes; conditions) |
| find-population-molecular-consequences | `find_population_molecular_consequences` | 1 |

That makes 14 sites. Each one currently computes the numerator from a list or set of patient IDs, hard-codes the denominator as `common.patients_db.count_documents({})`, and optionally emits `subject` parts.

## 2. Background (verified against code and the live DB)

- The app connects to MongoDB as a **readonly** user, so group data is loaded out of band with write credentials. This follows the same pattern as the patient and variant data.
- `Patients` had a duplicate document for `TCGA-DD-A1EH` (1116 documents, 1115 IDs). It was removed on 2026-11-03, so the denominator without `group` is now 1115. The expected test outputs still say 1116 (§5, step 0).
- Patient validation (`common.validate_subject`) returns **400**, while the README's status-code table says **404 Patient not found**. The IG asks for 404 for an unknown group (§4, decision 3).
- Integration tests run against the live Atlas DB and compare responses with JSON in `tests/expected_outputs/<operation>/`, so the groups used in tests have to exist in Atlas.
- The IG's example group `breast-cancer-trial-cohort` (40 members, including HG00403, HG00406 and HG00419) is illustrative only. It isn't loaded here, and HG00419 doesn't exist in our `Patients` collection.

## 3. Design

### 3.1 Group storage: `Group` collection in Mongo

`FHIRGenomicsData` has a `Group` collection, with one document per group:

| Card. | Field | Meaning |
| :--- | :--- | :--- |
| 1..1 | `_id` | Assigned by Mongo |
| 1..1 | `groupName` | Human-readable name of the group. Descriptive only, and not used by the operations. |
| 1..1 | `groupID` | The value supplied in the `group` input parameter. Enforced unique by an index. |
| 1..* | `memberID` | List of patient MRNs, i.e. `Patients.patientID` values |

```json
{
  "groupName": "HLA testing",
  "groupID": "hla-testing",
  "memberID": ["L2345", "NB6TK328", "XYZ123", "XYZ234", "XYZ345"]
}
```

- `memberID` has cardinality 1..*, so every group has at least one member and the empty-group case can't occur.
- Members don't have to exist in `Patients` (§4, decision 4).
- Indexes:
  - unique index on `groupID`: in place
  - no index on `memberID`. The operations never query by member. They fetch one group by `groupID` and filter in Python. An index on `memberID` would only help a reverse lookup ("which groups contain patient X?"), which nothing here needs.

**Loaded groups (2026-11-03).** These were pasted into Atlas via Insert Document, and verified with the readonly credential:

| groupID | groupName | Members |
| :--- | :--- | :--- |
| `hla-testing` | HLA testing | 5 |
| `pharmacogenomic-testing` | Pharmacogenomic testing | 1018 (incl. `PH1`…`PH1000`) |
| `germline-whole-exome-sequencing` | Germline Whole Exome Sequencing | 15 |
| `somatic-or-mixed-sequencing` | Somatic or Mixed Sequencing | 16 |

When the groups were prepared, a stray `patientID` entry (a header-row artifact) was dropped from `germline-whole-exome-sequencing`, and a duplicate `L2345` was removed from `somatic-or-mixed-sequencing`. All remaining members exist in `Patients`.

The source file `utilities/groups.json` was only a paste aid, and it is deleted before commit. Atlas is the source of truth.

Alternatives considered (not chosen):
- **Fetch a FHIR Group from an external FHIR server** (e.g. hapi.fhir.org, which the MTB app already uses for clinical data). This is closer to "a FHIR Group", but it adds a network dependency and a failure mode to every population query, and hapi's test data is volatile.
- **A JSON file in the repo.** This is simplest, but it's inconsistent with how all other server data lives in Mongo, and it needs a redeploy for every change.

### 3.2 Lookup helper (`app/common.py`)

```python
group_db = db.Group

def get_group_members(group):
    """Return the set of patientIDs in the group, or None if no group was supplied."""
    if not group:
        return None
    record = group_db.find_one({"groupID": group.strip()})
    if not record:
        abort(404, f"Group ({group}) not found.")
    return set(record["memberID"])
```

The helper is called once at the top of each population handler, before any query work, so an unknown group fails fast. Because it returns a set, a duplicate MRN within a group would count only once in the denominator.

### 3.3 Applying the filter: post-filter patient IDs

Each of the 14 sites already ends with a collection of patient IDs. The group filter intersects that collection with the group's `memberID` set before counting and listing.

This approach is the same at every site, it leaves the existing Mongo queries and the shared `common.query_*` helpers untouched, and it can't change results when `group` is omitted. The largest group has about 1,000 members, so the cost is negligible.

The alternative is to push `{"patientID": {"$in": memberID}}` into the Mongo `$match`. That would be more efficient at scale, but `query` is passed into roughly ten `common.query_*` helpers, and some of them restructure or `pop` keys from it. That makes it a riskier change across far more code.

### 3.4 Shared result helper (`app/common.py`)

The numerator/denominator/subject block that is copy-pasted at all 14 sites is replaced by one helper:

```python
def add_population_counts(parameter, patient_ids, group_members, includePatientList):
    if group_members is not None:
        patient_ids = set(patient_ids) & group_members
        denominator = len(group_members)
    else:
        denominator = patients_db.count_documents({})
    parameter["part"].append({"name": "numerator", "valueQuantity": {"value": len(patient_ids)}})
    parameter["part"].append({"name": "denominator", "valueQuantity": {"value": denominator}})
    if includePatientList:
        for patientID in sorted(patient_ids):
            parameter["part"].append({"name": "subject", "valueString": f"{patientID}"})
```

Each site then passes in its patient IDs: `[r["_id"] for r in query_results]`, `passed_patients`, `patient_ids`, and so on.

Without a group, the output must stay identical to today's. The refactor is done as its own step and checked against the existing expected-output tests before any group logic is added (§5, step 1).

### 3.5 OpenAPI spec ([app/api_spec.yml](app/api_spec.yml))

Add the following to each of the seven population operations, after `includePatientList`:

```yaml
- name: group
  in: query
  required: false
  schema:
    type: string
  description: Identifier of a FHIR Group known to the server. If supplied, the search is restricted to members of the group (numerator counts only group members, denominator is the number of members in the group, and the subject list contains only group members). If omitted, all patients in the repository are searched. An unknown group identifier results in a 404 error.
```

Also add a `404` response entry to these operations. The hosted API docs (redoc) are generated from this file, so they pick the change up automatically.

### 3.6 Unknown patient: 400 → 404

`common.validate_subject` ([app/common.py](app/common.py)) currently does `abort(400, f"Patient ({patient_id}) not found.")`. Change it to `abort(404, ...)`, so that an unknown patient matches the README (`404 | ERROR: Patient not found`) and is consistent with an unknown group.

- All 10 subject operations call this one function, so the code change is a single line, and no call site changes.
- Add a `404` response entry to the 10 subject operations in `app/api_spec.yml`.
- No existing test covers an unknown patient, so nothing breaks. New tests are listed in §3.8.
- The bundled MTB app doesn't depend on the 400. Any external client that does would see a behavior change, but the README already documents 404.

### 3.7 README

Add `404 | ERROR: Group not found` to the status-code table. The `404 | ERROR: Patient not found` row is already there and becomes accurate with §3.6.

### 3.8 Tests

New integration tests, using the loaded groups (§3.1):
- each population operation, run with `group=<loaded group>&includePatientList=true`, with new expected-output JSON files. For each operation, choose a group and query where some members match and some don't, so that both the numerator filtering and the denominator change are exercised.
- `group=does-not-exist` returns 404, tested on at least one genotype operation and one phenotype operation
- `subject=does-not-exist` returns 404, tested on at least one subject genotype operation and one subject phenotype operation
- all existing tests must pass unchanged

### 3.9 Out of scope

- Intensional groups (`Group.membership = 'definitional'`, R6). Groups are stored with an explicit member list only.
- Creating or editing groups through the API
- Using `group` in the MTB app

## 4. Decisions

All open questions are resolved.

1. **Storage:** a `Group` collection with `groupName` (1..1), `groupID` (1..1) and `memberID` (1..*) (§3.1).
2. **Identifier format:** `group` is the Group's business identifier (`Group.identifier`), not its logical id (`Group.id`). So the reference form `Group/xyz`, which points to a logical id, is **not** accepted. The value must exactly match `groupID` (after trimming whitespace). The `system|value` token form isn't supported either, because the schema has no identifier system. That could be added later if needed.
3. **Error status:** an unknown group returns 404 `Group not found`, per the IG. An unknown patient is changed from 400 to 404 `Patient not found` at the same time, for consistency and to match the README (§3.6).
4. **Denominator:** the denominator is the total count of `memberID`s, whether or not those members exist in `Patients`. Groups can be defined clinically, independently of the genomic repository. For example, the denominator might be "lung cancer patients" and the numerator "had EGFR testing", and not every group member will have genomic data. Members with no data in the repository simply never match, so they never appear in the numerator or the subject list.
5. **Refactor:** use the shared `add_population_counts` helper (§3.4). The refactor is a separate first step, with no group logic, verified against the existing tests before `group` is added.
6. **Groups:** four groups are loaded in Atlas (§3.1), with a unique index on `groupID` and no index on `memberID`.
7. **Empty group:** this can't happen, because `memberID` is 1..*.
8. **Duplicate members:** if an MRN appears twice in a group, it counts once (§3.2). None of the loaded groups has duplicates.

## 5. Implementation plan

Coding starts only on explicit go-ahead. Each step is reviewed before the next one begins.

0. **Make the expected-output tests able to fail.** `compare_actual_and_expected_output` in [tests/utilities.py](tests/utilities.py) asserts `json.dumps(sorted_expected) == json.dumps(sorted_expected)`, which compares the expected output with itself. So every expected-output test passes no matter what the server returns, unless the request itself errors. A strict run of the population tests (2026-11-03) gave **20 failed, 4 passed**: all 20 expected the old denominator of 1116 (it is now 1115, after the duplicate `TCGA-DD-A1EH` was removed from `Patients`), and `tx-implications` 3 and 5 and `dx-implications` 4 also gained newer matching patients (`L2345`, plus one more in dx 4). The fix:
   1. Compare `sorted_actual` with `sorted_expected` on that line.
   2. Run the full suite, subject operations included, and review every failure before regenerating anything.
   3. Regenerate the stale expected outputs with `OVERWRITE_TEST_EXPECTED_DATA=true python3 -m pytest`.
   4. Review `git diff tests/expected_outputs`. It should show only explained changes (the 1116→1115 denominators, the newer patients, and anything approved in step 0.2).

   This gives step 1 a test baseline that can actually catch a change in output.

   **Done (2026-11-03):** the full strict suite gave 23 failed, 74 passed. The failures were the 20 population tests above (the other new patient in dx 4 is `MRN1234`), plus 3 `find_study_metadata` tests (2, 3, 4), whose expected outputs had a stray leading space in the LOINC code `" LA26806-2"`. All 23 were regenerated. The diff contains only those changes, and all 97 tests pass.
1. **Refactor only (no behavior change).** Add `common.add_population_counts`, without the group branch active (`group_members=None` everywhere), and replace all 14 numerator/denominator/subject blocks with calls to it. Verify by running `python3 -m pytest`: every existing population test must pass with the expected outputs unchanged, and `git diff tests/expected_outputs` must be empty.
   **Done (2026-11-03):** `add_population_counts` was added to `app/common.py`, and all 14 blocks now call it (`app/endpoints.py`: 262 lines removed). Flake8 (CI settings) is clean. All 24 population tests pass, with no change to the expected outputs. Three blocks have no test coverage: the specific-variants per-item block, and the dx-implications variants and haplotypes paths. For those, the pre-refactor and refactored handlers were run side by side on the same queries, and their output was identical.
2. **Group lookup and unknown-patient status.** Add `common.group_db` and `common.get_group_members` (404 `Group not found`). Add the `group` parameter to the 7 handler signatures, and pass the members into `add_population_counts`. Change `common.validate_subject` from 400 to 404 (§3.6).
   **Done (2026-11-03):** `common.group_db`, `common.get_group_members` and the 404 in `validate_subject` were added. The 7 handlers take `group=None` and look the group up as their first statement. Until step 3, Connexion doesn't pass `group` over HTTP, so the handlers were called directly to check four operations against `germline-whole-exome-sequencing` and `pharmacogenomic-testing`: the denominator equals the group size, the numerator equals the number of subjects, and the subjects are group members and a subset of the unfiltered result. Unknown group → 404, unknown patient → 404 (over HTTP), and whitespace around `groupID` is trimmed. All 24 population tests still pass without `group`.
3. **Spec and README.** Add `group` and the 404 response to the 7 population operations, and the 404 response to the 10 subject operations, in `app/api_spec.yml`. Add the `404 | ERROR: Group not found` row to the README.
   **Done (2026-11-03):** `group` is in the 7 population operations, and a `404` response is in all 17 (7 population + 10 subject, including study metadata). The spec parses as YAML. Over HTTP, `$find-population-specific-haplotypes?haplotypes=CYP2C19` gives 1032/1115 without a group and 7/15 with `group=germline-whole-exome-sequencing`, and an unknown group gives 404. The README has the `Group not found` row. All 24 population tests pass.
4. **Tests.** Add the new integration tests (§3.8) and their expected outputs, using the loaded groups.
   **Done (2026-11-03):** 11 new tests were added. There are 7 group tests, one per population operation, with new expected outputs: specific-variants 4, structural-intersecting 4, structural-subsuming 3, specific-haplotypes 5, tx-implications 6, dx-implications 5 and molecular-consequences 2. Most use `somatic-or-mixed-sequencing`; haplotypes and dx use `germline-whole-exome-sequencing`. Queries were chosen where the group excludes some matching patients, where such a query exists. For subsuming and dx, no existing query/group combination does that, so those tests check only the denominator and the subject list. There are also 2 unknown-group 404 tests (specific-variants, tx-implications) and 2 unknown-patient 404 tests (find-subject-variants, find-subject-tx-implications). Each generated expected output was checked independently: the denominator equals the group size, the numerator equals the subject count, and the subjects are group members. The full suite passes: 108 tests.
5. **Documentation.**
   - **Wiki** (https://github.com/FHIR/genomics-operations/wiki, a separate git repo, `genomics-operations.wiki.git`):
     - *2. Getting Started*: describe the `Group` collection (`groupName`, `groupID`, `memberID`), list the available groups, and explain that a group's members don't all need genomic data in the repository (the denominator semantics).
     - *3. Examples and Exercises*: every section is named after an operation, so add `group` examples inside the relevant operation sections. Don't add a separate section.
     - *4. Software Utilities* and *6. Replicating the Reference Implementation*: no change expected. There's no loader script, and `mongodump`/`mongorestore` already copy every collection, including `Group`.
   - **Postman collection** (`FHIRGenomicsOperations.postman_collection.json`): add a `group` request for each population operation.
   - Wiki edits are drafted, then reviewed before being pushed.
   **Drafted (2026-11-03):** The Postman collection has one new `group` request per population folder (additions only; the file's formatting round-trips exactly). Wiki page 2 has a new "Patient Groups" section. On page 3, the group examples sit inside the `find-population-specific-haplotypes` section (CYP2C19 with pharmacogenomic-testing, unknown-group 404, and an exercise about the 5 members without CYP2C19 data) and the `find-population-tx-implications` section (NSCLC with somatic-or-mixed-sequencing). Pushed to `genomics-operations.wiki.git` (commits `e388844`, `43aa42b`, `b8b9e19`).
6. **Commit.** Delete `utilities/groups.json`, then commit the code, tests, spec, README, Postman collection and this design doc on `add-group-parameter`.
7. **After deployment: test the wiki examples.** The wiki examples call the production server (`fhir-gen-ops.herokuapp.com`). Until this branch is deployed there, production silently ignores `group` and returns unfiltered counts, so **don't test the wiki examples until this step**. Once deployed, run each new example on wiki pages 2 and 3 and confirm the results: NSCLC tx-implications with `somatic-or-mixed-sequencing` → 4 of 16; CYP2C19 with `pharmacogenomic-testing` → 1013 of 1018 (the exercise asks about the 5 missing members: HG00403, HG02657, NA19190, NA19239, NB6TK329, none of whom has a CYP2C19 genotype); unknown group → 404. Also run the new Postman requests.
