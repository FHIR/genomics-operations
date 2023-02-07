from flask import jsonify
from collections import OrderedDict
from .common import *


def find_subject_variants(
        subject, ranges, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, genomicSourceClass=None,
        includeVariants=False, includePhasing=False):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    ranges = list(map(get_range, ranges))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Genomics Build Presence
    genomics_build_presence = get_genomics_build_presence(query)

    # Chromosome To Ranges
    chromosome_to_ranges = get_chromosome_to_ranges(ranges, genomics_build_presence)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    # For each chromosome-range pair query the database for "variants" and
    # create a Parameters FHIR Object.
    for chrom in chromosome_to_ranges:
        parameter = OrderedDict()

        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "rangeItem",
          "valueString": f'{chrom["RefSeq"]}:{chrom["PGB"]["L"]}-{chrom["PGB"]["H"]}'
        })

        variant_q = []
        genomic_builds = [chrom["PGB"]["BUILD"]]

        query["$and"] = []
        query["$and"].append({"SVTYPE": { "$exists": False}})
        query["$and"].append({"$expr": {"$gte": [{"$subtract": [{"$add": [{"$strLenCP": "$REF"}, "$POS"]}, "$POS"]}, 1]}})
        query["$and"].append({"$or": [
            {
                "$and": [
                    {"POS": {"$lte": chrom["PGB"]["L"]}},
                    {"$expr": {"$gt":[{"$add": [{"$strLenCP": "$REF"}, "$POS"]}, chrom["PGB"]["L"]]}}
                ]
            },
            {
                "$and": [
                    {"POS": {"$gte": chrom["PGB"]["L"]}},
                    {"POS": {"$lt": chrom["PGB"]["H"]}}
                ]
            }
        ]})
        query["$and"].append({"CHROM": {"$eq": chrom["CHROM"]}})

        if chrom["OGB"] is not None:
            query["$and"][2]["$or"].extend(
                [{
                    "$and": [
                        {"POS": {"$lte": chrom["OGB"]["L"]}},
                        {"$expr": {"$gt":[{"$add": [{"$strLenCP": "$REF"}, "$POS"]}, chrom["OGB"]["L"]]}}
                    ]
                },
                {
                    "$and": [
                        {"POS": {"$gte": chrom["OGB"]["L"]}},
                        {"POS": {"$lt": chrom["OGB"]["H"]}}
                    ]
                }]
            )

            genomic_builds.append(chrom["OGB"]["BUILD"])

        query["genomicBuild"] = {"$in": genomic_builds}

        try:
            variant_q = variants_db.aggregate([{"$match": query}])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_variants query={query}")
            variant_q = []

        # Variants
        present = bool(variant_q)

        parameter["part"].append({
          "name": "presence",
          "valueBoolean": present
        })

        if present:
            if includeVariants:
                variant_fhir_profiles = []

                for record in variant_q:
                    ref_seq = get_ref_seq_by_chrom_and_build(record['genomicBuild'], record['CHROM'])
                    resource = create_fhir_variant_resource(record, ref_seq, subject)

                    variant_fhir_profiles.append(resource)

                if variant_fhir_profiles:
                    variant_fhir_profiles = sorted(variant_fhir_profiles, key=lambda d: d['id'])

                for resource in variant_fhir_profiles:
                    parameter["part"].append({
                      "name": "variant",
                      "resource": resource
                    })

                if includePhasing:
                    variantIDs = [str(v['_id']) for v in variant_q]
                    sequence_phase_profiles = []
                    sequence_phase_data = get_sequence_phase_data(subject)

                    for sq_data in sequence_phase_data:
                        if sq_data["variantID1"] in variantIDs and sq_data["variantID2"] in variantIDs:
                            sq_profile = create_sequence_phase_relationship(subject, sq_data)

                            sequence_phase_profiles.append(sq_profile)

                    if sequence_phase_profiles:
                        sequence_phase_profiles = sorted(sequence_phase_profiles, key=lambda d: d['id'])

                    for sq_profile in sequence_phase_profiles:
                        parameter["part"].append({
                          "name": "sequencePhaseRelationship",
                          "resource": sq_profile
                        })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_subject_specific_variants(
        subject, variants, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, genomicSourceClass=None):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    variants = list(map(get_variant, variants))

    # Query
    query = {}
    query["SVTYPE"] = {"$exists": False}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    genomics_build_presence = get_genomics_build_presence(query)

    for variant in variants:
        if not variant["GRCh37"] and genomics_build_presence["GRCh37"]:
            abort(422, f'Failed LiftOver. Variant: {variant["variant"]}')
        elif not variant["GRCh38"] and genomics_build_presence["GRCh38"]:
            abort(422, f'Failed LiftOver. Variant: {variant["variant"]}')

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    for variant in variants:
        parameter = OrderedDict()

        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "variantItem",
          "valueString": f"{variant['variant']}"
        })

        spdis = []

        if variant["GRCh37"]:
            spdis.append(variant["GRCh37"])
        if variant["GRCh38"]:
            spdis.append(variant["GRCh38"])

        query["SPDI"] = {"$in": spdis}

        try:
            variant_q = variants_db.aggregate([{"$match": query}])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_specific_variants query={query}")
            variant_q = []

        present = bool(variant_q)

        parameter["part"].append({
          "name": "presence",
          "valueBoolean": present
        })

        if present:
            variant_fhir_profiles = []

            for record in variant_q:
                ref_seq = get_ref_seq_by_chrom_and_build(record['genomicBuild'], record['CHROM'])
                resource = create_fhir_variant_resource(record, ref_seq, subject)

                variant_fhir_profiles.append(resource)

            if variant_fhir_profiles:
                variant_fhir_profiles = sorted(variant_fhir_profiles, key=lambda d: d['id'])

            for resource in variant_fhir_profiles:
                parameter["part"].append({
                  "name": "variant",
                  "resource": resource
                })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_subject_structural_intersecting_variants(
        subject, ranges, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, genomicSourceClass=None,
        includeVariants=False):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    ranges = list(map(get_range, ranges))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Genomics Build Presence
    genomics_build_presence = get_genomics_build_presence(query)

    # Chromosome To Ranges
    chromosome_to_ranges = get_chromosome_to_ranges(ranges, genomics_build_presence)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    # For each chromosome-range pair query the database for "variants" and
    # create a Parameters FHIR Object.
    for chrom in chromosome_to_ranges:
        parameter = OrderedDict()

        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "rangeItem",
          "valueString": f'{chrom["RefSeq"]}:{chrom["PGB"]["L"]}-{chrom["PGB"]["H"]}'
        })

        variant_q = []
        genomic_builds = [chrom["PGB"]["BUILD"]]

        query["$and"] = []
        query["$and"].append({"SVTYPE": { "$exists": True, "$ne": None}})
        query["$and"].append({"END": { "$exists": True, "$ne": None}})
        query["$and"].append({"$or": [
            {
                "$and": [
                    {"POS": {"$lte": chrom["PGB"]["L"]}},
                    {"END": {"$gte": chrom["PGB"]["L"]}}
                ]
            },
            {
                "$and": [
                    {"POS": {"$gte": chrom["PGB"]["L"]}},
                    {"POS": {"$lte": chrom["PGB"]["H"]}}
                ]
            }
        ]})
        query["$and"].append({"CHROM": {"$eq": chrom["CHROM"]}})

        if chrom["OGB"] is not None:
            query["$and"][2]["$or"].extend(
                [{
                    "$and": [
                        {"POS": {"$lte": chrom["OGB"]["L"]}},
                        {"END": {"$gt": chrom["OGB"]["L"]}}
                    ]
                },
                {
                    "$and": [
                        {"POS": {"$gte": chrom["OGB"]["L"]}},
                        {"POS": {"$lt": chrom["OGB"]["H"]}}
                    ]
                }]
            )

            genomic_builds.append(chrom["OGB"]["BUILD"])

        query["genomicBuild"] = {"$in": genomic_builds}

        try:
            variant_q = variants_db.aggregate([{"$match": query}])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_structural_intersecting_variants query={query}")
            variant_q = []

        # Variants
        present = bool(variant_q)

        parameter["part"].append({
          "name": "presence",
          "valueBoolean": present
        })

        if present:
            if includeVariants:
                variant_fhir_profiles = []

                for record in variant_q:
                    ref_seq = get_ref_seq_by_chrom_and_build(record['genomicBuild'], record['CHROM'])
                    resource = create_fhir_variant_resource(record, ref_seq, subject)

                    variant_fhir_profiles.append(resource)

                if variant_fhir_profiles:
                    variant_fhir_profiles = sorted(variant_fhir_profiles, key=lambda d: d['id'])

                for resource in variant_fhir_profiles:
                    parameter["part"].append({
                      "name": "variant",
                      "resource": resource
                    })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_subject_structural_subsuming_variants(
        subject, ranges, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, genomicSourceClass=None,
        includeVariants=False):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    ranges = list(map(get_range, ranges))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Genomics Build Presence
    genomics_build_presence = get_genomics_build_presence(query)

    # Chromosome To Ranges
    chromosome_to_ranges = get_chromosome_to_ranges(ranges, genomics_build_presence)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    # For each chromosome-range pair query the database for "variants" and
    # create a Parameters FHIR Object.
    for chrom in chromosome_to_ranges:
        parameter = OrderedDict()

        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "rangeItem",
          "valueString": f'{chrom["RefSeq"]}:{chrom["PGB"]["L"]}-{chrom["PGB"]["H"]}'
        })

        variant_q = []
        genomic_builds = [chrom["PGB"]["BUILD"]]

        query["$and"] = []
        query["$and"].append({"SVTYPE": { "$exists": True, "$ne": None}})
        query["$and"].append({"END": { "$exists": True, "$ne": None}})
        query["$and"].append({"$and": [
                {"POS": {"$lte": chrom["PGB"]["L"]}},
                {"END": {"$gte": chrom["PGB"]["H"]}}
            ]
        })
        query["$and"].append({"CHROM": {"$eq": chrom["CHROM"]}})

        if chrom["OGB"] is not None:
            query["$and"][2]["$and"].extend(
                [
                    {"POS": {"$lte": chrom["OGB"]["L"]}},
                    {"END": {"$gte": chrom["OGB"]["H"]}}
                ]
            )

            genomic_builds.append(chrom["OGB"]["BUILD"])

        query["genomicBuild"] = {"$in": genomic_builds}

        try:
            variant_q = variants_db.aggregate([{"$match": query}])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_structural_subsuming_variants query={query}")
            variant_q = []

        # Variants
        present = bool(variant_q)

        parameter["part"].append({
          "name": "presence",
          "valueBoolean": present
        })

        if present:
            if includeVariants:
                variant_fhir_profiles = []

                for record in variant_q:
                    ref_seq = get_ref_seq_by_chrom_and_build(record['genomicBuild'], record['CHROM'])
                    resource = create_fhir_variant_resource(record, ref_seq, subject)

                    variant_fhir_profiles.append(resource)

                if variant_fhir_profiles:
                    variant_fhir_profiles = sorted(variant_fhir_profiles, key=lambda d: d['id'])

                for resource in variant_fhir_profiles:
                    parameter["part"].append({
                      "name": "variant",
                      "resource": resource
                    })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_subject_haplotypes(
        subject, genes, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, genomicSourceClass=None):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    genes = list(map(get_gene, genes))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    # if genomicSourceClass:
    #     genomicSourceClass = genomicSourceClass.strip().lower()
    #     query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    for gene in genes:
        parameter = OrderedDict()

        parameter["name"] = "haplotypes"
        parameter["part"] = []
        parameter["part"].append({
          "name": "geneItem",
          "valueString": f"{gene['gene']}"
        })

        query['$or'] = []

        if gene['isSystem']:
            query['$or'].append({'geneCode': {"$eq": gene['gene']}})
        else:
            query['$or'].append({'$or': [
                    {'geneCode': {'$regex': ".*"+str(gene['gene']).replace('*', r'\*')+".*"}},
                    {'geneDesc': {'$regex': ".*"+str(gene['gene']).replace('*', r'\*')+".*"}}
                ]})

        try:
            gene_q = genotypes_db.aggregate([{"$match": query}])
            gene_q = list(gene_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_haplotypes query={query}")
            gene_q = []

        genotype_profiles = []
        for qresult in gene_q:
            # haplotype_profile = create_haplotype_profile(qresult, subject, "")
            genotype_profile = create_genotype_profile(qresult, subject, [])

            genotype_profiles.append(genotype_profile)

            # parameter["part"].append({
            # "name": "haplotype",
            # "resource": haplotype_profile
            # })

        if genotype_profiles:
            genotype_profiles = sorted(genotype_profiles, key=lambda d: d['id'])

        for genotype_profile in genotype_profiles:
            parameter["part"].append({
            "name": "genotype",
            "resource": genotype_profile
            })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_subject_specific_haplotypes(
        subject, haplotypes, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, genomicSourceClass=None):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    haplotypes = list(map(get_haplotype, haplotypes))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    # if genomicSourceClass:
    #     genomicSourceClass = genomicSourceClass.strip().lower()
    #     query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    for haplotype in haplotypes:
        parameter = OrderedDict()

        parameter["name"] = "haplotypes"
        parameter["part"] = []
        parameter["part"].append({
          "name": "haplotypeItem",
          "valueString": f"{haplotype['haplotype']}"
        })

        if haplotype['isSystem']:
            query['$and'] = [
                {'genotypeCode': {"$eq": haplotype['haplotype']}},
                {'genotypeCodeSystem': {"$eq": haplotype['system']}}
            ]
        else:
            query['$or'] = [
                    {'genotypeCode': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}},
                    {'genotypeDesc': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}}
                ]

        try:
            haplotype_q = genotypes_db.aggregate([{"$match": query}])
            haplotype_q = list(haplotype_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_specific_haplotypes query={query}")
            haplotype_q = []

        present = bool(haplotype_q)

        parameter["part"].append({
          "name": "presence",
          "valueBoolean": present
        })

        if present:
            genotype_profiles = []
            for qresult in haplotype_q:
                # haplotype_profile = create_haplotype_profile(qresult, subject, "")
                genotype_profile = create_genotype_profile(qresult, subject, [])

                genotype_profiles.append(genotype_profile)

                # parameter["part"].append({
                # "name": "haplotype",
                # "resource": haplotype_profile
                # })
                
            if genotype_profiles:
                genotype_profiles = sorted(genotype_profiles, key=lambda d: d['id'])

            for genotype_profile in genotype_profiles:
                parameter["part"].append({
                "name": "genotype",
                "resource": genotype_profile
                })

            result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_subject_tx_implications(
        subject, variants=None, ranges=None, haplotypes=None, treatments=None, conditions=None,
        testIdentifiers=None, testDateRange=None, specimenIdentifiers=None,
        genomicSourceClass=None):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    if variants and ranges:
        abort(400, f"You cannot supply both 'variants' and 'ranges'.")

    if (variants or ranges) and haplotypes:
        abort(400, f"You cannot supply both ('variants' or 'ranges') and 'haplotypes'.")
        
    if not variants and not conditions and not treatments and not haplotypes and not ranges:
        return jsonify({"resourceType":"Parameters"})
    
    normalized_haplotype_list = []
    if haplotypes:
        normalized_haplotype_list = list(map(get_haplotype, haplotypes))    
    
    condition_code_list = []
    if conditions:
        condition_code_list = list(map(get_condition, conditions))
    
    treatment_code_list = []
    if treatments:
        treatment_code_list = list(map(get_treatment, treatments))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    normalized_variants = []
    if ranges:
        ranges = list(map(get_range, ranges))
        get_lift_over_range(ranges)
        variants = get_variants(ranges, query)
        if not variants:
            return jsonify({"resourceType":"Parameters"})
        normalized_variants = [{variant["BUILD"]: variant["SPDI"]} for variant in variants]

    if variants and not ranges:
        normalized_variants = list(map(get_variant, variants))

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    if normalized_variants:
        if not ranges:
            genomics_build_presence = get_genomics_build_presence(query)

            for normalizedVariant in normalized_variants:
                if not normalizedVariant["GRCh37"] and genomics_build_presence["GRCh37"]:
                    abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')
                elif not normalizedVariant["GRCh38"] and genomics_build_presence["GRCh38"]:
                    abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')

        query_results = query_CIVIC_by_variants(normalized_variants, condition_code_list, treatment_code_list, query)

        for res in query_results:
            if res["txImplicationMatches"]:
                ref_seq = get_ref_seq_by_chrom_and_build(res['genomicBuild'], res['CHROM'])
            for implication in res["txImplicationMatches"]:
                parameter = OrderedDict()
                parameter["name"] = "implications"
                parameter["part"] = []

                implication_profile = create_tx_implication_profile_civic(implication, subject, [str(res['_id'])])
                parameter["part"].append({
                "name": "implication",
                "resource": implication_profile
                })
                resource = create_fhir_variant_resource(res, ref_seq, subject)

                add_variation_id(resource, implication["variationID"])

                parameter["part"].append({
                "name": "variant",
                "resource": resource
                })

                result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        else:
            result["parameter"] = sorted(result["parameter"], key=lambda d: d['part'][0]['resource']['id'])

        return jsonify(result)
    
    if haplotypes:
        if genomicSourceClass:
            query.pop("genomicSourceClass")
        
        query_results = query_PharmGKB_by_haplotypes(normalized_haplotype_list, treatment_code_list, query)
        print(query_results)
        for res in query_results:
            for implication in res["txImplicationMatches"]:
                parameter = OrderedDict()
                parameter["name"] = "implications"
                parameter["part"] = []

                implication_profile = create_tx_implication_profile_pharmgkb(implication, subject, [str(res['_id'])])
                parameter["part"].append({
                "name": "implication",
                "resource": implication_profile
                })
                
                # haplotype_profile = create_haplotype_profile(res, subject, res["UUID"])

                # parameter["part"].append({
                # "name": "haplotype",
                # "resource": haplotype_profile
                # })
                
                genotype_profile = create_genotype_profile(res, subject, [str(res['_id'])])

                add_variation_id(genotype_profile, implication["variationID"])

                parameter["part"].append({
                "name": "genotype",
                "resource": genotype_profile
                })

                result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        else:
            result["parameter"] = sorted(result["parameter"], key=lambda d: d['part'][0]['resource']['id'])

        return jsonify(result)
    
    if treatments:
        query_results_PGKB = query_PharmGKB_by_treatments(condition_code_list, treatment_code_list, query)
        query_results_CIViC = query_CIVIC_by_condition(condition_code_list, treatment_code_list, query)
        
        for res in query_results_PGKB:
            parameter = OrderedDict()
            parameter["name"] = "implications"
            parameter["part"] = []

            implication_profile = create_tx_implication_profile_pharmgkb(res, subject, [str(i['_id']) for i in res["patientMatches"]])
            parameter["part"].append({
            "name": "implication",
            "resource": implication_profile
            })
            genotype_profiles = []
            for genItem in res["patientMatches"]:
                
                # haplotype_profile = create_haplotype_profile(genItem, subject, genItem["UUID"])

                # parameter["part"].append({
                # "name": "haplotype",
                # "resource": haplotype_profile
                # })
                
                genotype_profile = create_genotype_profile(genItem, subject, [str(genItem['_id'])])

                add_variation_id(genotype_profile, res["variationID"])

                genotype_profiles.append(genotype_profile)

            if genotype_profiles:
                genotype_profiles = sorted(genotype_profiles, key=lambda d: d['id'])

            for genotype_profile in genotype_profiles:
                parameter["part"].append({
                "name": "genotype",
                "resource": genotype_profile
                })

            result["parameter"].append(parameter)
        
        for res in query_results_CIViC:
            parameter = OrderedDict()
            parameter["name"] = "implications"
            parameter["part"] = []

            implication_profile = create_tx_implication_profile_civic(res, subject, [str(i['_id']) for i in res["patientMatches"]])
            parameter["part"].append({
            "name": "implication",
            "resource": implication_profile
            })

            variant_fhir_profiles = []
            for varItem in res["patientMatches"]:
                ref_seq = get_ref_seq_by_chrom_and_build(varItem['genomicBuild'], varItem['CHROM'])
                resource = create_fhir_variant_resource(varItem, ref_seq, subject)

                add_variation_id(resource, res["variationID"])

                variant_fhir_profiles.append(resource)

            if variant_fhir_profiles:
                variant_fhir_profiles = sorted(variant_fhir_profiles, key=lambda d: d['id'])

            for resource in variant_fhir_profiles:
                parameter["part"].append({
                "name": "variant",
                "resource": resource
                })

            result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        else:
            result["parameter"] = sorted(result["parameter"], key=lambda d: d['part'][0]['resource']['id'])
            
        return jsonify(result)
    
    if conditions:
        query_results = query_CIVIC_by_condition(condition_code_list, treatment_code_list, query)
        for res in query_results:
            parameter = OrderedDict()
            parameter["name"] = "implications"
            parameter["part"] = []

            implication_profile = create_tx_implication_profile_civic(res, subject, [str(i['_id']) for i in res["patientMatches"]])
            parameter["part"].append({
            "name": "implication",
            "resource": implication_profile
            })

            variant_fhir_profiles = []
            for varItem in res["patientMatches"]:
                ref_seq = get_ref_seq_by_chrom_and_build(varItem['genomicBuild'], varItem['CHROM'])
                resource = create_fhir_variant_resource(varItem, ref_seq, subject)

                add_variation_id(resource, res["variationID"])

                variant_fhir_profiles.append(resource)
                
            if variant_fhir_profiles:
                variant_fhir_profiles = sorted(variant_fhir_profiles, key=lambda d: d['id'])

            for resource in variant_fhir_profiles:
                parameter["part"].append({
                "name": "variant",
                "resource": resource
                })

            result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        else:
            result["parameter"] = sorted(result["parameter"], key=lambda d: d['part'][0]['resource']['id'])

        return jsonify(result)


def find_subject_dx_implications(
        subject, variants=None, ranges=None, conditions=None, testIdentifiers=None,
        testDateRange=None, specimenIdentifiers=None, genomicSourceClass=None):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    if not variants and not conditions and not ranges:
        return jsonify({"resourceType":"Parameters"})

    if variants and ranges:
        abort(400, f"You cannot supply both 'variants' and 'ranges'.")

    condition_code_list = []
    if conditions:
        condition_code_list = list(map(get_condition, conditions))

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    normalized_variants = []
    if ranges:
        ranges = list(map(get_range, ranges))
        get_lift_over_range(ranges)
        variants = get_variants(ranges, query)
        if not variants:
            return jsonify({"resourceType":"Parameters"})
        normalized_variants = [{variant["BUILD"]: variant["SPDI"]} for variant in variants]
    
    if variants and not ranges:
        normalized_variants = list(map(get_variant, variants))

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    if normalized_variants:
        if not ranges:
            genomics_build_presence = get_genomics_build_presence(query)

            for normalizedVariant in normalized_variants:
                if not normalizedVariant["GRCh37"] and genomics_build_presence["GRCh37"]:
                    abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')
                elif not normalizedVariant["GRCh38"] and genomics_build_presence["GRCh38"]:
                    abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')

        query_results = query_clinvar_by_variants(normalized_variants, condition_code_list, query)

        for res in query_results:
            if res["dxImplicationMatches"]:
                ref_seq = get_ref_seq_by_chrom_and_build(res['genomicBuild'], res['CHROM'])
            for implication in res["dxImplicationMatches"]:
                parameter = OrderedDict()
                parameter["name"] = "implications"
                parameter["part"] = []

                implication_profile = create_dx_implication_profile(implication, subject, [str(res['_id'])])

                parameter["part"].append({
                "name": "implication",
                "resource": implication_profile
                })
                resource = create_fhir_variant_resource(res, ref_seq, subject)

                add_variation_id(resource, implication["variationID"])

                parameter["part"].append({
                "name": "variant",
                "resource": resource
                })

                result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        else:
            result["parameter"] = sorted(result["parameter"], key=lambda d: d['part'][0]['resource']['id'])

        return jsonify(result)
    
    if conditions:
        query_results = query_clinvar_by_condition(condition_code_list, query)

        for res in query_results:
            parameter = OrderedDict()
            parameter["name"] = "implications"
            parameter["part"] = []

            implication_profile = create_dx_implication_profile(res, subject, [str(i['_id']) for i in res["patientMatches"]])
            parameter["part"].append({
            "name": "implication",
            "resource": implication_profile
            })

            for varItem in res["patientMatches"]:
                ref_seq = get_ref_seq_by_chrom_and_build(varItem['genomicBuild'], varItem['CHROM'])
                resource = create_fhir_variant_resource(varItem, ref_seq, subject)

                add_variation_id(resource, res["variationID"])

                parameter["part"].append({
                "name": "variant",
                "resource": resource
                })

            result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        else:
            result["parameter"] = sorted(result["parameter"], key=lambda d: d['part'][0]['resource']['id'])

        return jsonify(result)


def find_study_metadata(
        subject, testIdentifiers=None, testDateRange=None,
        specimenIdentifiers=None, ranges=None):

    # Parameters
    subject = subject.strip()
    validate_subject(subject)

    if ranges:
        ranges = list(map(get_range, ranges))
        ranges = merge_ranges(ranges)

    # Query
    query = {}

    # date query
    if testDateRange:
        testDateRange = list(map(get_date, testDateRange))
        query["testDate"] = {}

        for date_range in testDateRange:
            query["testDate"][date_range['OPERATOR']] = date_range['DATE']

    # Subject Query
    query["patientID"] = {"$eq": subject}

    # testIdentifiers query
    if testIdentifiers:
        testIdentifiers = list(map(str.strip, testIdentifiers))
        query["testID"] = {"$in": testIdentifiers}

    # specimenIdentifiers query
    if specimenIdentifiers:
        specimenIdentifiers = list(map(str.strip, specimenIdentifiers))
        query["specimenID"] = {"$in": specimenIdentifiers}

    # Genomics Build Presence
    genomics_build_presence = get_genomics_build_presence_tests_db(query)

    # Chromosome To Ranges
    chromosome_to_ranges = get_chromosome_to_ranges(ranges, genomics_build_presence)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    try:
        identified_tests = tests_db.aggregate([{"$match": query}])
        identified_tests = list(identified_tests)
    except Exception as e:
        print(f"DEBUG: Error{e} under find_study_metadata query={query}")
        identified_tests = []

    for identified_test in identified_tests:
        parameter = OrderedDict()

        parameter["name"] = "tests"
        parameter["part"] = []

        parameter["part"].append({
          "name": "testID",
          "valueString": f"{identified_test['testID']}"
        })

        parameter["part"].append({
          "name": "testDate",
          "valueString": f"{identified_test['testDate']}"
        })

        parameter["part"].append({
          "name": "specimenId",
          "valueString": f"{identified_test['specimenID']}"
        })

        parameter["part"].append({
          "name": "genomicBuild",
          "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                               "code": f"{GENOMIC_BUILD_TO_CODE[identified_test['genomicBuild']]}",
                                               "display": f"{identified_test['genomicBuild']}"}]}
        })

        if 'dnaChangeType' in identified_test:
            parameter["part"].append({"name": "dnaChangeType",
                                      "valueCodeableConcept": {"coding": []}
                })
            for dct in identified_test["dnaChangeType"]:
                if dct.strip() in DNA_CHANGE_TYPE_TO_CODE:
                    parameter["part"][-1]["valueCodeableConcept"]["coding"].append(
                        {"system": "http://sequenceontology.org",
                         "code": f"{DNA_CHANGE_TYPE_TO_CODE[dct.strip()]}",
                         "display": f"{dct.strip()}"}
                     )

        if chromosome_to_ranges:
            intersected_region_studied = []
            intersected_uncallable_regions = []

            region_studied_present = 'studiedRegion' in identified_test
            uncallable_region_present = 'uncallableRegion' in identified_test
            for chromosome_to_range in chromosome_to_ranges:
                if region_studied_present:
                    if chromosome_to_range["PGB"]["BUILD"] == identified_test['genomicBuild']:
                        get_intersected_regions(identified_test["studiedRegion"],
                            chromosome_to_range["PGB"]["BUILD"],
                            chromosome_to_range["CHROM"],
                            chromosome_to_range["PGB"]["L"],
                            chromosome_to_range["PGB"]["H"],
                            intersected_region_studied)
                    elif chromosome_to_range["OGB"] and chromosome_to_range["OGB"]["BUILD"] == identified_test['genomicBuild']:
                        get_intersected_regions(identified_test["studiedRegion"],
                            chromosome_to_range["OGB"]["BUILD"],
                            chromosome_to_range["CHROM"],
                            chromosome_to_range["OGB"]["L"],
                            chromosome_to_range["OGB"]["H"],
                            intersected_region_studied)
                    else:
                        abort(422, f'Failed LiftOver ({chromosome_to_range["PGB"]["RefSeq"]}:{chromosome_to_range["PGB"]["L"]}-{chromosome_to_range["PGB"]["L"]})')

                if uncallable_region_present:
                    if chromosome_to_range["PGB"]["BUILD"] == identified_test['genomicBuild']:
                        get_intersected_regions(identified_test["uncallableRegion"],
                            chromosome_to_range["PGB"]["BUILD"],
                            chromosome_to_range["CHROM"],
                            chromosome_to_range["PGB"]["L"],
                            chromosome_to_range["PGB"]["H"],
                            intersected_uncallable_regions)
                    elif chromosome_to_range["OGB"] and chromosome_to_range["OGB"]["BUILD"] == identified_test['genomicBuild']:
                        get_intersected_regions(identified_test["uncallableRegion"],
                            chromosome_to_range["OGB"]["BUILD"],
                            chromosome_to_range["CHROM"],
                            chromosome_to_range["OGB"]["L"],
                            chromosome_to_range["OGB"]["H"],
                            intersected_uncallable_regions)
                    else:
                        abort(422, f'Failed LiftOver ({chromosome_to_range["PGB"]["RefSeq"]}:{chromosome_to_range["PGB"]["L"]}-{chromosome_to_range["PGB"]["L"]})')

            if region_studied_present:
                parameter["part"].append({
                          "name": "regionStudied",
                          "valueString": f"{intersected_region_studied}"
                        })
            else:
                parameter["part"].append({
                  "name": "regionStudied",
                  "valueString": "unknown"
                })

            if uncallable_region_present:
                parameter["part"].append({
                      "name": "uncallableRegions",
                      "valueString": f"{intersected_uncallable_regions}"
                    })
            else:
                parameter["part"].append({
                  "name": "uncallableRegions",
                  "valueString": "unknown"
                })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_population_specific_variants(
    variants, genomicSourceClass=None, includePatientList=None):

    # Parameters
    variants = list(map(lambda x: x.strip().split(","), variants))
    for i in range(len(variants)):
        variants[i] = list(map(get_variant, variants[i]))

    # Query
    query = {}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    variantItem = []
    normalized_variant_lists = []
    for normalized_variant_list in variants:
        variantList = []
        providedVariantList = []

        for normalizedVariant in normalized_variant_list:
            if normalizedVariant["GRCh37"]:
                variantList.append(normalizedVariant["GRCh37"])
            if normalizedVariant["GRCh38"]:
                variantList.append(normalizedVariant["GRCh38"])

            providedVariantList.append(normalizedVariant["variant"])

        variantItem.append('|'.join(providedVariantList))
        normalized_variant_lists.append(variantList)

    variantItem = ' AND '.join(variantItem)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    if len(variants) > 1:
        parameter = OrderedDict()

        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "variantItem",
          "valueString": f"{variantItem}"
        })

        all_patients = []
        for varList in normalized_variant_lists:
            query["SPDI"] = {"$in": varList}

            try:
                variant_q = variants_db.aggregate([
                    {"$match": query},
                    {'$group': {'_id': '$patientID'}}
                ])
                variant_q = list(variant_q)
            except Exception as e:
                print(f"DEBUG: Error{e} under find_population_specific_variants query={query}")
                variant_q = []

            patients = []

            for patientID in variant_q:
                patients.append(patientID['_id'])

            all_patients.append(set(patients))

        passed_patients = set.intersection(*all_patients)

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {"value": len(passed_patients)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            for patientID in sorted(passed_patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f"{patientID}"
                })

        result["parameter"].append(parameter)

    else:
        for varItem in variants[0]:
            parameter = OrderedDict()

            parameter["name"] = "variants"
            parameter["part"] = []
            parameter["part"].append({
              "name": "variantItem",
              "valueString": f"{varItem['variant']}"
            })

            varList = []

            if 'GRCh37' in varItem:
                varList.append(varItem["GRCh37"])

            if 'GRCh38' in varItem:
                varList.append(varItem["GRCh38"])

            query["SPDI"] = {"$in": varList}
            try:
                variant_q = variants_db.aggregate([
                    {"$match": query},
                    {'$group': {'_id': '$patientID'}}
                ])
                variant_q = list(variant_q)
            except Exception as e:
                print(f"DEBUG: Error{e} under find_population_specific_variants query={query}")
                variant_q = []

            parameter["part"].append({
              "name": "numerator",
              "valueQuantity": {"value": len(variant_q)}
            })

            parameter["part"].append({
              "name": "denominator",
              "valueQuantity": {"value": patients_db.count_documents({})}
            })

            patients = []

            for patientID in variant_q:
                patients.append(f"{patientID['_id']}")

            if includePatientList:
                for patientID in sorted(patients):
                    parameter["part"].append({
                      "name": "subject",
                      "valueString": f"{patientID}"
                    })

            result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_population_structural_intersecting_variants(
    ranges, genomicSourceClass=None, includePatientList=None):

    # Parameters
    ranges = list(map(get_range, ranges))

    # Query
    query = {}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    chromosome_to_ranges = get_chromosome_to_ranges(ranges, None)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    for chrom in chromosome_to_ranges:
        parameter = OrderedDict()
        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "rangeItem",
          "valueString": f'{chrom["RefSeq"]}:{chrom["PGB"]["L"]}-{chrom["PGB"]["H"]}'
        })

        query["$and"] = []
        query["$and"].append({"SVTYPE": { "$exists": True, "$ne": None}})
        query["$and"].append({"END": { "$exists": True, "$ne": None}})
        query["$and"].append({"$or": [
            {
                "$and": [
                    {"POS": {"$lte": chrom["PGB"]["L"]}},
                    {"END": {"$gte": chrom["PGB"]["L"]}}
                ]
            },
            {
                "$and": [
                    {"POS": {"$gte": chrom["PGB"]["L"]}},
                    {"POS": {"$lte": chrom["PGB"]["H"]}}
                ]
            }
        ]})
        query["$and"].append({"CHROM": {"$eq": chrom["CHROM"]}})

        if chrom["OGB"] is not None:
            query["$and"][2]["$or"].extend(
                [{
                    "$and": [
                        {"POS": {"$lte": chrom["OGB"]["L"]}},
                        {"END": {"$gte": chrom["OGB"]["L"]}}
                    ]
                },
                {
                    "$and": [
                        {"POS": {"$gte": chrom["OGB"]["L"]}},
                        {"POS": {"$lte": chrom["OGB"]["H"]}}
                    ]
                }]
            )

        try:
            variant_q = variants_db.aggregate([
                {"$match": query},
                {'$group': {'_id': '$patientID'}}
            ])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_population_structural_intersecting_variants query={query}")
            variant_q = []

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(variant_q)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            patients = []
            for patientID in variant_q:
                patients.append(f'{patientID["_id"]}')

            for patientID in sorted(patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_population_structural_subsuming_variants(
    ranges, genomicSourceClass=None, includePatientList=None):

    # Parameters
    ranges = list(map(get_range, ranges))

    # Query
    query = {}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    chromosome_to_ranges = get_chromosome_to_ranges(ranges, None)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    for chrom in chromosome_to_ranges:
        parameter = OrderedDict()
        parameter["name"] = "variants"
        parameter["part"] = []
        parameter["part"].append({
          "name": "rangeItem",
          "valueString": f'{chrom["RefSeq"]}:{chrom["PGB"]["L"]}-{chrom["PGB"]["H"]}'
        })

        query["$and"] = []
        query["$and"].append({"SVTYPE": { "$exists": True, "$ne": None}})
        query["$and"].append({"END": { "$exists": True, "$ne": None}})
        query["$and"].append({"$and": [
                {"POS": {"$lte": chrom["PGB"]["L"]}},
                {"END": {"$gte": chrom["PGB"]["H"]}}
            ]
        })
        query["$and"].append({"CHROM": {"$eq": chrom["CHROM"]}})

        if chrom["OGB"] is not None:
            query["$and"][2]["$and"].extend(
                [
                    {"POS": {"$lte": chrom["OGB"]["L"]}},
                    {"END": {"$gte": chrom["OGB"]["H"]}}
                ]
            )

        try:
            variant_q = variants_db.aggregate([
                {"$match": query},
                {'$group': {'_id': '$patientID'}}
            ])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_population_structural_subsuming_variants query={query}")
            variant_q = []

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(variant_q)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            patients = []
            for patientID in variant_q:
                patients.append(f'{patientID["_id"]}')

            for patientID in sorted(patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_population_specific_haplotypes(
    haplotypes, genomicSourceClass=None, includePatientList=None):

    # Parameters
    haplotypes = list(map(lambda x: x.strip().split(","), haplotypes))
    for i in range(len(haplotypes)):
        haplotypes[i] = list(map(get_haplotype, haplotypes[i]))

    # Query
    query = {}

    # Genomic Source Class Query
    # if genomicSourceClass:
    #     genomicSourceClass = genomicSourceClass.strip().lower()
    #     query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    haplotypeItem = []
    normalizedHaplotypesLists = []
    for normalized_haplotype_list in haplotypes:
        haplotypeList = []
        providedHaplotypeList = []

        for normalizedHaplotype in normalized_haplotype_list:
            haplotypeList.append(normalizedHaplotype)
            providedHaplotypeList.append(normalizedHaplotype['haplotype'])

        haplotypeItem.append('|'.join(providedHaplotypeList))
        normalizedHaplotypesLists.append(haplotypeList)

    haplotypeItem = ' AND '.join(haplotypeItem)

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    if len(haplotypes) > 1:

        parameter = OrderedDict()

        parameter["name"] = "haplotypes"
        parameter["part"] = []
        parameter["part"].append({
          "name": "haplotypeItem",
          "valueString": f"{haplotypeItem}"
        })

        all_patients = []
        for hapList in normalizedHaplotypesLists:

            for haplotype in hapList:
                if haplotype['isSystem']:
                    query['$and'] = [
                        {'genotypeCode': {"$eq": haplotype['haplotype']}},
                        {'genotypeCodeSystem': {"$eq": haplotype['system']}}
                    ]
                else:
                    query['$or'] = [
                        {'genotypeCode': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}},
                        {'genotypeDesc': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}}
                    ]

            try:
                haplotype_q = genotypes_db.aggregate([
                    {"$match": query},
                    {'$group': {'_id': '$patientID'}}
                ])
                haplotype_q = list(haplotype_q)
            except Exception as e:
                print(f"DEBUG: Error{e} under find_population_specific_haplotypes query={query}")
                haplotype_q = []

            patients = []

            for patientID in haplotype_q:
                patients.append(patientID['_id'])

            all_patients.append(set(patients))

        passed_patients = set.intersection(*all_patients)

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {"value": len(passed_patients)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            for patientID in sorted(passed_patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f"{patientID}"
                })

        result["parameter"].append(parameter)

    else:
        for hapItem in haplotypes[0]:
            parameter = OrderedDict()

            parameter["name"] = "variants"
            parameter["part"] = []
            parameter["part"].append({
              "name": "variantItem",
              "valueString": f"{hapItem['haplotype']}"
            })

            if hapItem['isSystem']:
                query['$and'] = [
                    {'genotypeCode': {"$eq": hapItem['haplotype']}},
                    {'genotypeCodeSystem': {"$eq": hapItem['system']}}
                ]
            else:
                query['$or'] = [
                        {'genotypeCode': {'$regex': ".*"+str(hapItem['haplotype']).replace('*', r'\*')+".*"}},
                        {'genotypeDesc': {'$regex': ".*"+str(hapItem['haplotype']).replace('*', r'\*')+".*"}}
                    ]

            try:
                haplotype_q = genotypes_db.aggregate([
                    {"$match": query},
                    {'$group': {'_id': '$patientID'}}
                ])
                haplotype_q = list(haplotype_q)
            except Exception as e:
                print(f"DEBUG: Error{e} under find_population_specific_haplotypes query={query}")
                haplotype_q = []

            parameter["part"].append({
              "name": "numerator",
              "valueQuantity": {"value": len(haplotype_q)}
            })

            parameter["part"].append({
              "name": "denominator",
              "valueQuantity": {"value": patients_db.count_documents({})}
            })

            patients = []

            for patientID in haplotype_q:
                patients.append(f"{patientID['_id']}")

            if includePatientList:
                for patientID in sorted(patients):
                    parameter["part"].append({
                      "name": "subject",
                      "valueString": f"{patientID}"
                    })

            result["parameter"].append(parameter)

    if not result["parameter"]:
        result.pop("parameter")

    return jsonify(result)


def find_population_tx_implications(
        variants=None, haplotypes=None, treatments=None, conditions=None,
        genomicSourceClass=None, includePatientList=None):

    # Parameters
    if not variants and not conditions and not treatments and not haplotypes:
        return jsonify({"resourceType":"Parameters"})

    if variants and haplotypes:
        return jsonify({"resourceType":"Parameters"})

    if variants:
        variants = list(map(get_variant, variants))
    
    condition_code_list = []
    if conditions:
        condition_code_list = list(map(get_condition, conditions))

    normalized_haplotype_list = []
    if haplotypes:
        normalized_haplotype_list = list(map(get_haplotype, haplotypes))

    treatment_code_list = []
    if treatments:
        treatment_code_list = list(map(get_treatment, treatments))

    # Query
    query = {}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    if variants:
        genomics_build_presence = get_genomics_build_presence(query)
        for normalizedVariant in variants:
            if not normalizedVariant["GRCh37"] and genomics_build_presence["GRCh37"]:
                abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')
            elif not normalizedVariant["GRCh38"] and genomics_build_presence["GRCh38"]:
                abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')

        query_results = query_CIVIC_by_variants(variants, condition_code_list, treatment_code_list, query, True)

        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(query_results)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            patients = []
            for patientID in query_results:
                patients.append(f'{patientID["_id"]}')

            for patientID in sorted(patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")

        return jsonify(result)
    
    if haplotypes:
        if genomicSourceClass:
            query.pop("genomicSourceClass")

        query_results = query_PharmGKB_by_haplotypes(normalized_haplotype_list, [], query, True)
        
        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(query_results)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            patients = []
            for patientID in query_results:
                patients.append(f'{patientID["_id"]}')

            for patientID in sorted(patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")

        return jsonify(result)

    if treatments:
        query_results_PGKB = query_PharmGKB_by_treatments(condition_code_list, treatment_code_list, query)
        query_results_CIViC = query_CIVIC_by_condition(condition_code_list, treatment_code_list, query)
        
        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        patient_ids = []

        if query_results_PGKB:
            for query_result in query_results_PGKB:
                for patientID in query_result['patientMatches']:
                    patient_ids.append(patientID["patientID"])

        if query_results_CIViC:
            for query_result in query_results_CIViC:
                for patientID in query_result['patientMatches']:
                    patient_ids.append(patientID["patientID"])

        patient_ids = list(set(patient_ids))

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(patient_ids)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            for patientID in sorted(patient_ids):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })
        
        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
            
        return jsonify(result)

    if conditions:
        query_results = query_CIVIC_by_condition(condition_code_list, treatment_code_list, query)

        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        patient_ids = []

        if query_results:
            for query_result in query_results:
                for patientID in query_result['patientMatches']:
                    patient_ids.append(patientID["patientID"])

        patient_ids = list(set(patient_ids))

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(patient_ids)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            for patientID in sorted(patient_ids):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")

        return jsonify(result)


def find_population_dx_implications(
        variants=None, haplotypes=None, conditions=None, genomicSourceClass=None,
        includePatientList=None):

    # Parameters
    if not variants and not conditions and not haplotypes:
        return jsonify({"resourceType":"Parameters"})

    if variants and haplotypes:
        return jsonify({"resourceType":"Parameters"})

    if variants:
        variants = list(map(get_variant, variants))

    condition_code_list = []
    if conditions:
        condition_code_list = list(map(get_condition, conditions))

    normalized_haplotype_list = []
    if haplotypes:
        normalized_haplotype_list = list(map(get_haplotype, haplotypes))

    # Query
    query = {}

    # Genomic Source Class Query
    if genomicSourceClass:
        genomicSourceClass = genomicSourceClass.strip().lower()
        query["genomicSourceClass"] = {"$eq": genomicSourceClass}

    # Result Object
    result = OrderedDict()
    result["resourceType"] = "Parameters"
    result["parameter"] = []

    if variants:
        genomics_build_presence = get_genomics_build_presence(query)

        for normalizedVariant in variants:
            if not normalizedVariant["GRCh37"] and genomics_build_presence["GRCh37"]:
                abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')
            elif not normalizedVariant["GRCh38"] and genomics_build_presence["GRCh38"]:
                abort(422, f'Failed LiftOver. Variant: {normalizedVariant["variant"]}')

        query_results = query_clinvar_by_variants(variants, condition_code_list, query, True)

        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(query_results)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            patients = []
            for patientID in query_results:
                patients.append(f'{patientID["_id"]}')

            for patientID in sorted(patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        return jsonify(result)
        
    if haplotypes:
        if genomicSourceClass:
            query.pop("genomicSourceClass")

        query_results = query_PharmGKB_by_haplotypes(normalized_haplotype_list, [], query, True)
        
        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(query_results)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            patients = []
            for patientID in query_results:
                patients.append(f'{patientID["_id"]}')

            for patientID in sorted(patients):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")
        return jsonify(result)

    if conditions:
        query_results = query_clinvar_by_condition(condition_code_list, query)

        parameter = OrderedDict()
        parameter["name"] = "implications"
        parameter["part"] = []

        patient_ids = []

        for query_result in query_results:
            for patientID in query_result['patientMatches']:
                patient_ids.append(patientID["patientID"])

        patient_ids = list(set(patient_ids))

        parameter["part"].append({
          "name": "numerator",
          "valueQuantity": {'value': len(patient_ids)}
        })

        parameter["part"].append({
          "name": "denominator",
          "valueQuantity": {"value": patients_db.count_documents({})}
        })

        if includePatientList:
            for patientID in sorted(patient_ids):
                parameter["part"].append({
                  "name": "subject",
                  "valueString": f'{patientID}'
                })

        result["parameter"].append(parameter)

        if not result["parameter"]:
            result.pop("parameter")

        return jsonify(result)
