import { FhirResponse, Component, FhirObservation, Coding } from "@/utils/fhirInterfaces";

type ComponentMap = Record<string, string>;

type ExtraFieldsExtractor<T> = (resource: FhirObservation) => Partial<T>;

export const processFhirResponse = <T>(
    response: FhirResponse,
    componentMap: ComponentMap,
    variantIdSystem?: string,
    extraFieldsExtractor?: ExtraFieldsExtractor<T>
): T[] => {
    if (!response?.parameter) return [];

    return response.parameter
        .filter(param => {
            if (!param?.resource?.meta?.profile?.[0]) return false;
            const profile = param.resource.meta.profile[0];
            return (
                (param.name === 'implication' && (
                    profile.includes('diagnostic-implication') ||
                    profile.includes('therapeutic-implication')
                )) ||
                (param.name === 'consequence' && profile.includes('molecular-consequence'))
            );
        })
        .map(param => {
            const { component = [] } = param.resource;

            const componentFields = Object.fromEntries(
                Object.entries(componentMap).map(([key, code]) => [
                    key,
                    findComponentValue(component, code)
                ])
            );

            const extraFields = extraFieldsExtractor ? extraFieldsExtractor(param.resource) : {};

            const result = {
                ...componentFields,
                ...extraFields,
            };

            return result as T;
        });
};

const componentsWithCode = (components: Component[], code: string) =>
    components.filter(c => c.code.coding?.some(coding => coding.code === code));

const joinedDisplays = (components: Component[], code: string, separator: string) =>
    componentsWithCode(components, code)
        .map(c => c.valueCodeableConcept?.coding?.[0]?.display || c.valueCodeableConcept?.text || '')
        .filter(Boolean)
        .join(separator);

/**
 * Every coding of every component with the code, e.g. all of an implication's cancer types (one component each)
 */
export const findComponentCodings = (components: Component[] | undefined, code: string): Coding[] =>
    componentsWithCode(components ?? [], code).flatMap(c => c.valueCodeableConcept?.coding ?? []);

export const findComponentValue = (components: Component[], code: string): string => {
    // A combination therapy has one medication-assessed component per drug: show them all, joined with +
    if (code === '51963-7') {
        return joinedDisplays(components, code, ' + ');
    }

    // An implication can apply to several cancer types, one component each: show them all, joined with ;
    if (code === '81259-4') {
        return joinedDisplays(components, code, '; ');
    }

    const component = components.find(c =>
        c.code.coding?.some(coding => coding.code === code)
    );

    if (code === '53037-8' || code === 'feature-consequence') {
        return component?.valueCodeableConcept?.coding?.[0]?.display || '';
    }

    // HGVS and the transcript come as text (the c.HGVS was read from valueString before, so it was always empty)
    if (code === '48004-6' || code === '51958-7' || code === '48005-3' || code === '93044-6') {
        return component?.valueCodeableConcept?.text || '';
    }

    // For predicted-therapeutic-implication, check valueCodeableConcept.text first, then fallback to valueString
    if (code === 'predicted-therapeutic-implication') {
        return component?.valueCodeableConcept?.text || component?.valueString || '';
    }

    return component?.valueString || '';
};
