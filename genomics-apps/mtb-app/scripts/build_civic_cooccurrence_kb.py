"""Build public/data/CIViC_Cooccurrence_KB.csv from the CIViC GraphQL API.

The file lists CIViC molecular profiles that describe co-occurring variants, one row per
predictive evidence item. The app flags a profile when every component is among the CIViC
identifiers returned with the patient's therapeutic implications: either the component's CIViC
variant ID or its single-variant molecular profile ID (Cat-VRS results identify the latter).

Run by hand from the mtb-app folder, then commit the CSV and review the git diff:
    python scripts/build_civic_cooccurrence_kb.py

Uses only the Python standard library.
"""

import csv
import json
import sys
import time
import urllib.request
from datetime import date
from pathlib import Path

CIVIC_GRAPHQL_URL = "https://civicdb.org/api/graphql"
OUTPUT_PATH = Path(__file__).resolve().parent.parent / "public" / "data" / "CIViC_Cooccurrence_KB.csv"

# Filters
# - Profiles combine two or more variants with AND only. Profiles containing OR are skipped,
#   because "every component present" would wrongly require every alternative.
# - Profiles with a NOT / Wildtype part are skipped (no clean definition of wildtype yet).
# - Every component must be a gene variant; fusions and non-gene factors are skipped.
# - Deprecated profiles are skipped.
# - Only predictive evidence with status ACCEPTED or SUBMITTED is included.
INCLUDED_EVIDENCE_STATUSES = {"ACCEPTED", "SUBMITTED"}
INCLUDED_EVIDENCE_TYPES = {"PREDICTIVE"}
PAGE_SIZE = 100
DETAIL_BATCH_SIZE = 25

CSV_COLUMNS = [
    "mp_id", "mp_name", "component_variant_ids", "component_profile_ids", "genes", "evidence_id", "disease", "therapies",
    "interaction_type", "level", "direction", "significance", "status", "civic_url", "retrieved_date",
]

PROFILE_LIST_QUERY = """
query($after: String) {
  molecularProfiles(first: %d, after: $after) {
    pageInfo { hasNextPage endCursor }
    nodes { id isMultiVariant deprecated }
  }
}
""" % PAGE_SIZE

PROFILE_DETAIL_QUERY = """
query($ids: [Int!]) {
  molecularProfiles(ids: $ids, first: %d) {
    nodes {
      id
      name
      deprecated
      parsedName { __typename ... on MolecularProfileTextSegment { text } }
      variants { __typename id name feature { name } ... on GeneVariant { singleVariantMolecularProfile { id } } }
      evidenceItems(first: 200) {
        nodes {
          id
          status
          evidenceType
          evidenceLevel
          evidenceDirection
          significance
          therapyInteractionType
          disease { name }
          therapies { name }
        }
      }
    }
  }
}
""" % DETAIL_BATCH_SIZE


def run_query(query, variables):
    payload = json.dumps({"query": query, "variables": variables}).encode()
    request = urllib.request.Request(CIVIC_GRAPHQL_URL, payload, {"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                body = json.load(response)
            if body.get("errors"):
                raise RuntimeError(body["errors"])
            return body["data"]
        except Exception as error:  # network hiccups: retry with backoff
            if attempt == 3:
                raise
            print(f"  retrying after error: {error}", file=sys.stderr)
            time.sleep(2 ** attempt)


def list_multi_variant_profile_ids():
    ids, after = [], None
    while True:
        page = run_query(PROFILE_LIST_QUERY, {"after": after})["molecularProfiles"]
        ids += [int(node["id"]) for node in page["nodes"] if node["isMultiVariant"] and not node["deprecated"]]
        if not page["pageInfo"]["hasNextPage"]:
            return ids
        after = page["pageInfo"]["endCursor"]


def skip_reason(profile):
    operators = {segment["text"].strip().upper() for segment in profile["parsedName"]
                 if segment["__typename"] == "MolecularProfileTextSegment"}
    variant_names = {variant["name"].lower() for variant in profile["variants"]}
    if profile["deprecated"]:
        return "deprecated"
    if "AND" not in operators:
        return "no AND"
    if "OR" in operators:
        return "contains OR"
    if "NOT" in operators or "wildtype" in variant_names:
        return "NOT / wildtype"
    if any(variant["__typename"] != "GeneVariant" or "::" in variant["feature"]["name"] for variant in profile["variants"]):
        return "fusion or non-gene component"
    if len({variant["id"] for variant in profile["variants"]}) < 2:
        return "fewer than 2 variants"
    return None


def main():
    retrieved = date.today().isoformat()
    print("Listing CIViC molecular profiles...")
    profile_ids = list_multi_variant_profile_ids()
    print(f"  {len(profile_ids)} multi-variant profiles")

    rows, skipped, kept_profiles = [], {}, 0
    for start in range(0, len(profile_ids), DETAIL_BATCH_SIZE):
        batch = profile_ids[start:start + DETAIL_BATCH_SIZE]
        for profile in run_query(PROFILE_DETAIL_QUERY, {"ids": batch})["molecularProfiles"]["nodes"]:
            reason = skip_reason(profile)
            evidence = [item for item in profile["evidenceItems"]["nodes"]
                        if item["evidenceType"] in INCLUDED_EVIDENCE_TYPES and item["status"] in INCLUDED_EVIDENCE_STATUSES]
            if not reason and not evidence:
                reason = "no predictive evidence"
            if reason:
                skipped[reason] = skipped.get(reason, 0) + 1
                continue
            kept_profiles += 1
            components = sorted({(int(variant["id"]), int(variant["singleVariantMolecularProfile"]["id"])) for variant in profile["variants"]})
            genes = sorted({variant["feature"]["name"] for variant in profile["variants"]})
            for item in evidence:
                rows.append({
                    "mp_id": profile["id"],
                    "mp_name": profile["name"],
                    "component_variant_ids": "; ".join(str(variant_id) for variant_id, _ in components),
                    # Same order as component_variant_ids
                    "component_profile_ids": "; ".join(str(profile_id) for _, profile_id in components),
                    "genes": "; ".join(genes),
                    "evidence_id": item["id"],
                    "disease": (item["disease"] or {}).get("name", ""),
                    # Separated by "; "; the app joins them according to interaction_type
                    "therapies": "; ".join(therapy["name"] for therapy in item["therapies"]),
                    "interaction_type": item["therapyInteractionType"] or "",
                    "level": item["evidenceLevel"] or "",
                    "direction": item["evidenceDirection"] or "",
                    "significance": item["significance"] or "",
                    "status": item["status"],
                    "civic_url": f"https://civicdb.org/molecular-profiles/{profile['id']}/summary",
                    "retrieved_date": retrieved,
                })

    rows.sort(key=lambda row: (int(row["mp_id"]), int(row["evidence_id"])))
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} evidence rows for {kept_profiles} profiles to {OUTPUT_PATH}")
    print("Skipped profiles: " + ", ".join(f"{reason}: {count}" for reason, count in sorted(skipped.items())))


if __name__ == "__main__":
    main()
