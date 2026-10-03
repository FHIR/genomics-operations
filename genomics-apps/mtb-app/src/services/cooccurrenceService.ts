// CO-OCCURRING VARIANTS (CIViC molecular profiles that combine variants with AND)
import Papa from "papaparse";
import { Variant } from '@/types/variants';

const COOCCURRENCE_KB_PATH = '/data/CIViC_Cooccurrence_KB.csv';

type CooccurrenceKbRow = {
    mp_id: string;
    mp_name: string;
    component_variant_ids: string;
    component_profile_ids: string;
    genes: string;
    evidence_id: string;
    disease: string;
    therapies: string;
    interaction_type: string;
    level: string;
    direction: string;
    significance: string;
    civic_url: string;
};

export interface CooccurrenceEvidence {
    evidenceId: string;
    disease: string;
    therapies: string;
    interactionType: string;
    level: string;
    direction: string;
    significance: string;
}

export interface CooccurrenceComponent {
    variantId: string;
    profileId?: string; // The component's single-variant molecular profile
}

export interface CooccurrenceProfile {
    id: string;
    name: string;
    components: CooccurrenceComponent[];
    genes: string[];
    url: string;
    evidence: CooccurrenceEvidence[];
}

export interface CooccurrenceMatch {
    profile: CooccurrenceProfile;
    variantIds: string[]; // Patient variants (Variant.id) that supply the profile's components
    componentVariantIds: string[][]; // For each component, the patient variants that supply it
}

const splitList = (value?: string) => (value ?? '').split(';').map((item) => item.trim()).filter(Boolean);

let profilesPromise: Promise<CooccurrenceProfile[]> | null = null;

// Loads the knowledge base once per page view
export function loadCooccurrenceProfiles() {
    if (!profilesPromise) {
        profilesPromise = new Promise<CooccurrenceProfile[]>((resolve, reject) => {
            Papa.parse<CooccurrenceKbRow>(COOCCURRENCE_KB_PATH, {
                header: true,
                download: true,
                skipEmptyLines: true,
                complete: (result) => {
                    const profiles = new Map<string, CooccurrenceProfile>();

                    result.data.forEach((row) => {
                        const id = row.mp_id?.trim();
                        if (!id) {
                            return;
                        }

                        const profileIds = splitList(row.component_profile_ids);
                        const profile = profiles.get(id) ?? {
                            id,
                            name: row.mp_name.trim(),
                            components: splitList(row.component_variant_ids).map((variantId, index) => ({ variantId, profileId: profileIds[index] })),
                            genes: splitList(row.genes),
                            url: row.civic_url.trim(),
                            evidence: [],
                        };

                        profile.evidence.push({
                            evidenceId: row.evidence_id.trim(),
                            disease: row.disease.trim(),
                            therapies: row.therapies.trim(),
                            interactionType: row.interaction_type.trim(),
                            level: row.level.trim(),
                            direction: row.direction.trim(),
                            significance: row.significance.trim(),
                        });
                        profiles.set(id, profile);
                    });

                    resolve(Array.from(profiles.values()));
                },
                error: (error) => {
                    profilesPromise = null;
                    reject(error);
                },
            });
        });
    }

    return profilesPromise;
}

// A profile matches when every component is among the CIViC identifiers returned with the
// patient's therapeutic implications: the component's CIViC variant ID, or its single-variant
// molecular profile ID (Cat-VRS results are identified by molecular profile)
export function findCooccurrenceMatches(variants: Variant[], profiles: CooccurrenceProfile[]): CooccurrenceMatch[] {
    const variantIdsByCivicKey = new Map<string, Set<string>>();
    const addKey = (key: string, variantId: string) => {
        const variantIds = variantIdsByCivicKey.get(key) ?? new Set<string>();
        variantIds.add(variantId);
        variantIdsByCivicKey.set(key, variantIds);
    };

    variants.forEach((variant) => {
        if (!variant.id) {
            return;
        }

        (variant.txImplications ?? []).forEach((implication) => {
            (implication.civicVariantIds ?? []).forEach((civicId) => addKey(`variant:${civicId}`, variant.id as string));
            (implication.civicProfileIds ?? []).forEach((civicId) => addKey(`profile:${civicId}`, variant.id as string));
        });
    });

    const variantIdsForComponent = (component: CooccurrenceComponent) => [
        ...(variantIdsByCivicKey.get(`variant:${component.variantId}`) ?? []),
        ...(component.profileId ? variantIdsByCivicKey.get(`profile:${component.profileId}`) ?? [] : []),
    ];

    return profiles
        .filter((profile) => profile.components.every((component) => variantIdsForComponent(component).length > 0))
        .map((profile) => {
            const componentVariantIds = profile.components.map((component) => [...new Set(variantIdsForComponent(component))]);
            return { profile, componentVariantIds, variantIds: [...new Set(componentVariantIds.flat())] };
        });
}

// Matches still complete when only the selected variants are considered: every component needs at
// least one selected variant. Each returned match lists only the selected variants that supply it.
export function restrictMatchesToVariants(matches: CooccurrenceMatch[], variantIds: Set<string>): CooccurrenceMatch[] {
    return matches.flatMap((match) => {
        const componentVariantIds = match.componentVariantIds.map((ids) => ids.filter((id) => variantIds.has(id)));
        if (componentVariantIds.some((ids) => ids.length === 0)) {
            return [];
        }
        return [{ ...match, componentVariantIds, variantIds: [...new Set(componentVariantIds.flat())] }];
    });
}

export function groupMatchesByVariant(matches: CooccurrenceMatch[]) {
    const byVariant: Record<string, CooccurrenceProfile[]> = {};

    matches.forEach((match) => {
        match.variantIds.forEach((variantId) => {
            byVariant[variantId] = [...(byVariant[variantId] ?? []), match.profile];
        });
    });

    return byVariant;
}

const levelRank = (level: string) => {
    const rank = 'ABCDE'.indexOf(level.trim().toUpperCase());
    return rank === -1 ? 5 : rank;
};

// Evidence listed by level, A (strongest) through E; same level by CIViC evidence ID
export function sortEvidenceByLevel(evidence: CooccurrenceEvidence[]) {
    return [...evidence].sort((left, right) => levelRank(left.level) - levelRank(right.level)
        || Number(left.evidenceId) - Number(right.evidenceId));
}

// Same wording as the app's therapeutic implications:
// "{subject} {implication} to {medication} in {phenotype} [{evidence level}]"
const SIGNIFICANCE_TEXT: Record<string, string> = {
    SENSITIVITYRESPONSE: 'Sensitivity/Response',
    RESISTANCE: 'Resistance',
    REDUCED_SENSITIVITY: 'Reduced Sensitivity',
    ADVERSE_RESPONSE: 'Adverse Response',
};

const LEVEL_TEXT: Record<string, string> = {
    A: 'A (Validated association)',
    B: 'B (Clinical evidence)',
    C: 'C (Case study)',
    D: 'D (Preclinical evidence)',
    E: 'E (Inferential association)',
};

// COMBINATION: "A + B"; SUBSTITUTES: "A or B" / "A, B or C"; SEQUENTIAL: "A, then B"
export function formatTherapies(therapies: string, interactionType: string) {
    const names = splitList(therapies);

    if (names.length < 2) {
        return names[0] ?? '';
    }

    if (interactionType === 'SUBSTITUTES') {
        return `${names.slice(0, -1).join(', ')} or ${names[names.length - 1]}`;
    }

    if (interactionType === 'SEQUENTIAL') {
        return names.join(', then ');
    }

    return names.join(' + ');
}

export function describeCooccurrenceEvidence(evidence: CooccurrenceEvidence) {
    const verb = evidence.direction === 'DOES_NOT_SUPPORT' ? 'Do Not Support' : 'Support';
    const significance = SIGNIFICANCE_TEXT[evidence.significance] ?? evidence.significance;
    const parts = ['Co-occurring Variants', verb, significance];

    const therapies = formatTherapies(evidence.therapies, evidence.interactionType);

    if (therapies) {
        parts.push('to', therapies);
    }

    if (evidence.disease) {
        parts.push('in', evidence.disease);
    }

    if (evidence.level) {
        parts.push(`[${LEVEL_TEXT[evidence.level] ?? evidence.level}]`);
    }

    return parts.join(' ');
}
