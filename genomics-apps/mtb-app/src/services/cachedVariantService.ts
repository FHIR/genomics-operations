import { Variant } from '@/types/variants';
import { findSubjectVariants } from './variantService';
import { processDxImplications } from './dxService';
import { processTxImplications } from './txService';
import { processMolecularConsequences } from './mcService';
import cacheService from './cacheService';

interface ProcessedVariant extends Variant {
  id: string;
  searchId: number;
}

/**
 * Cached version of findSubjectVariants that includes all implications
 * This function will check the cache first and return cached results if available
 */
export const findSubjectVariantsWithCache = async (
  range: string,
  subjectId: string,
  searchId: number
): Promise<ProcessedVariant[]> => {
  // Check cache first
  const cachedEntry = await cacheService.getCachedResults(range, subjectId);
  if (cachedEntry) {
    console.log(`Using cached data for ${range}`);

    // Convert cached data back to ProcessedVariant format
    return cachedEntry.response.processedVariants.map(cachedVariant => ({
      ...cachedVariant,
      id: `${range}-${cachedVariant.range}-${cachedVariant.variant}-${searchId}`,
      searchId,
      isLoading: false
    }));
  }
  
  console.log(` Fetching new data for ${range}`);
  
  // If not in cache, fetch from API
  const variants = await findSubjectVariants(range, subjectId);
  
  if (variants.length === 0) {
    // Cache empty results too to avoid repeated API calls
    await cacheService.setCachedResults(range, subjectId, []);
    return [];
  }
  
  // Process implications for all variants
  const processedVariants: ProcessedVariant[] = [];
  
  for (const variant of variants) {
    try {
      // Fetch all implications concurrently
      const [dx, tx, mc] = await Promise.allSettled([
        processDxImplications(variant.variant, subjectId),
        processTxImplications(variant.variant, subjectId),
        processMolecularConsequences(variant.variant, subjectId)
      ]);

      // Extract results and handle failures gracefully
      const dxImplications = dx.status === 'fulfilled' ? dx.value : [];
      const txImplications = tx.status === 'fulfilled' ? tx.value : [];
      const molecularConsequences = mc.status === 'fulfilled' ? mc.value : [];

      // Log individual failures
      if (dx.status === 'rejected') {
        console.error(`DX implications failed for ${variant.variant}:`, dx.reason);
      }
      if (tx.status === 'rejected') {
        console.error(`TX implications failed for ${variant.variant}:`, tx.reason);
      }
      if (mc.status === 'rejected') {
        console.error(`Molecular consequences failed for ${variant.variant}:`, mc.reason);
      }

      // Create processed variant
      const processedVariant: ProcessedVariant = {
        ...variant,
        id: `${range}-${variant.range}-${variant.variant}-${searchId}`,
        dxImplications,
        txImplications,
        molecularConsequences,
        searchId,
        isLoading: false
      };

      processedVariants.push(processedVariant);

    } catch (error) {
      console.error(`Error processing variant ${variant.variant}:`, error);

      // Create variant with error state
      const errorVariant: ProcessedVariant = {
        ...variant,
        id: `${range}-${variant.range}-${variant.variant}-${searchId}`,
        dxImplications: [],
        txImplications: [],
        molecularConsequences: [],
        searchId,
        isLoading: false,
        error: error instanceof Error ? error.message : 'Unknown error'
      };

      processedVariants.push(errorVariant);
    }
  }

  // Convert processed variants to cache format
  const cachedVariants = processedVariants.map(pv => ({
    id: pv.id,
    range: pv.range,
    variant: pv.variant,
    dxImplications: pv.dxImplications,
    txImplications: pv.txImplications,
    molecularConsequences: pv.molecularConsequences
  }));

  // Cache the results for future use
  await cacheService.setCachedResults(range, subjectId, cachedVariants);
  
  console.log(`Cached results for ${range}`);
  
  return processedVariants;
};
