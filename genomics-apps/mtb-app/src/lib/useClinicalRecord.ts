import { useCallback, useEffect, useState } from 'react';
import { ClinicalRecord, fetchClinicalRecord } from '@/services/clinicalService';

export type ClinicalState =
  | { status: 'idle' }
  | { status: 'loading' }
  | { status: 'none' }
  | { status: 'error'; message: string }
  | { status: 'loaded'; record: ClinicalRecord };

// Wait for typing in the custom MRN box to pause before querying the FHIR server
const DEBOUNCE_MS = 400;

/**
 * Loads the clinical record for an MRN. Responses for an MRN that is no longer selected are ignored.
 */
export function useClinicalRecord(mrn: string) {
  const [state, setState] = useState<ClinicalState>({ status: 'idle' });
  const [retryCount, setRetryCount] = useState(0);

  useEffect(() => {
    const trimmedMrn = mrn.trim();
    if (!trimmedMrn) {
      setState({ status: 'idle' });
      return;
    }

    let isCurrent = true;
    setState({ status: 'loading' });
    const timer = window.setTimeout(() => {
      fetchClinicalRecord(trimmedMrn)
        .then((record) => {
          if (isCurrent) {
            setState(record ? { status: 'loaded', record } : { status: 'none' });
          }
        })
        .catch((error) => {
          console.error(`Error loading clinical data for ${trimmedMrn}:`, error);
          if (isCurrent) {
            setState({ status: 'error', message: error instanceof Error ? error.message : String(error) });
          }
        });
    }, DEBOUNCE_MS);

    return () => {
      isCurrent = false;
      window.clearTimeout(timer);
    };
  }, [mrn, retryCount]);

  const retry = useCallback(() => setRetryCount((count) => count + 1), []);

  return { state, retry };
}
