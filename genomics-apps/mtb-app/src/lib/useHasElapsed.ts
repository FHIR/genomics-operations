import { useEffect, useState } from 'react';

/**
 * Returns true once `thresholdMs` has passed since `startedAt`.
 * Pass undefined for `startedAt` when nothing is being timed.
 */
export function useHasElapsed(startedAt: number | undefined, thresholdMs: number) {
    const [hasElapsed, setHasElapsed] = useState(false);

    useEffect(() => {
        if (startedAt === undefined) {
            setHasElapsed(false);
            return;
        }

        const remainingMs = thresholdMs - (Date.now() - startedAt);

        if (remainingMs <= 0) {
            setHasElapsed(true);
            return;
        }

        setHasElapsed(false);
        const timeoutId = setTimeout(() => setHasElapsed(true), remainingMs);
        return () => clearTimeout(timeoutId);
    }, [startedAt, thresholdMs]);

    return hasElapsed;
}
