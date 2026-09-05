import { FhirResponse, Component, FhirObservation } from "@/utils/fhirInterfaces";

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

export const findComponentValue = (components: Component[], code: string): string => {
    const component = components.find(c =>
        c.code.coding?.some(coding => coding.code === code)
    );

    if (code === '53037-8' || code === '81259-4' || code === 'feature-consequence' || code === '51963-7') {
        return component?.valueCodeableConcept?.coding?.[0]?.display || '';
    }

    if (code === '48005-3' || code === '93044-6') {
        return component?.valueCodeableConcept?.text || '';
    }

    // For predicted-therapeutic-implication, check valueCodeableConcept.text first, then fallback to valueString
    if (code === 'predicted-therapeutic-implication') {
        return component?.valueCodeableConcept?.text || component?.valueString || '';
    }

    return component?.valueString || '';
};
