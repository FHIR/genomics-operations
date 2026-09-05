export interface FhirResponse {
    parameter: FhirParameter[];
}

export interface FhirParameter {
    name: string;
    resource: FhirObservation;
}
export interface FhirObservation {
    resourceType: string;
    id?: string;
    component: Component[];
    interpretation?: {
        text?: string;
        coding?: {
            system?: string;
            code?: string;
            display?: string;
        }[];
    }[];
    identifier?: { system: string; value: string }[];
    subject?: { reference: string };
    derivedFrom?: { reference: string }[];
    meta?: { profile?: string[] };
    valueCodeableConcept?: CodeableConcept;
    valueString?: string;
}

export interface Component {
    code: CodeableConcept;
    valueCodeableConcept?: CodeableConcept;
    valueString?: string;
    valueRange?: {
        low?: {
            value: number;
        };
        high?: {
            value: number;
        };
    };
}

export interface CodeableConcept {
    coding: Coding[];
    text?: string;
}

export interface Coding {
    system?: string;
    code: string;
    display?: string;
}
