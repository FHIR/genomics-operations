# HL7 FHIR Genomics Operations - Reference Implementation
Source code for a [public reference implementation](https://fhir-gen-ops.herokuapp.com/) of [HL7 FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html).

Please refer to [project Wiki page](https://github.com/FHIR/genomics-operations/wiki) for details of this reference implementation, including how to replicate.

For additional information on the operations and the reference implementation, please see 

hi there!!

cer patient to available clinical trials based on identified somatic variants; screening for actionable hereditary conditions; identifying a risk for adverse medication reactions based on pharmacogenomic variants; updating a patient's risk as knowledge of their variants evolves; and more. A goal for FHIR Genomics operations is to ultimately support any and all of these clinical scenarios.

## Scope

In scope are clinical genomics operations. In the future, operations supporting variant calling and annotation, and knowledge base lookups may be added. We further categorize clinical genomics operations along two orthogonal axes - subject vs. population, and genotype vs. phenotype. For example, the '**find-subject-variants**' operation retrieves genotype information for a single subject; whereas the '**find-population-tx-implications**' retrieves a count or list of patients having specific phenotypes (such as being intermediate metabolizers of clopidogrel).

## Response

All operations return a JSON output. However, if an invalid request is submitted, or some other error occurs, a JSON response is returned in the following format:

```javascript
{
  "type": string,
  "title": string,
  "detail": string,
  "status": int
}
```

## Status Codes

The operations return the following status codes:

| Status Code | Description |
| :--- | :--- |
| 200 | `Successfully executed request` |
| 400 | `ERROR: Invalid query parameters` |
| 404 | `ERROR: Patient not found` |
| 422 | `ERROR: Failed LiftOver` |
| 500 | `INTERNAL SERVER ERROR` |

## Testing

To run the [integration tests](https://github.com/FHIR/genomics-operations/tree/main/tests), you can use the VS Code Testing functionality which should discover them automatically. You can also
run `python3 -m pytest` from the terminal to execute them all.

Additionally, since the tests run against the Mongo DB database, if you need to update the test data in this repo, you
can run `OVERWRITE_TEST_EXPECTED_DATA=true python3 -m pytest` from the terminal and then create a pull request with the
changes.
