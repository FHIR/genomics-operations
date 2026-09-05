import { DxImplication } from './dxService';
import { ProcessedTxImplication } from './txService';
import { MolecularConsequence } from './mcService';



interface ProcessedVariantCache {
  id: string;
  range: string;
  variant: string;
  dxImplications: DxImplication[];
  txImplications: ProcessedTxImplication[];
  molecularConsequences: MolecularConsequence[];
}

interface CacheEntry {
  timestamp: string;
  date: string; // YYYY-MM-DD format for quick date filtering
  query: {
    range: string;
    subjectId: string;
  };
  response: {
    processedVariants: ProcessedVariantCache[];
  };
}

interface CacheLog {
  entries: CacheEntry[];
}

class CacheService {
  private static instance: CacheService;
  private cache: Map<string, CacheEntry> = new Map();
  private isInitialized = false;

  private constructor() {}

  static getInstance(): CacheService {
    if (!CacheService.instance) {
      CacheService.instance = new CacheService();
    }
    return CacheService.instance;
  }

  /**
   * Generate cache key from query parameters
   */
  private getCacheKey(range: string, subjectId = 'L2345'): string {
    return `${range}:${subjectId}`;
  }

  /**
   * Check if cache entry is from today
   */
  private isFromToday(entry: CacheEntry): boolean {
    const today = new Date().toISOString().split('T')[0]; // YYYY-MM-DD
    return entry.date === today;
  }

  /**
   * Initialize cache by loading from localStorage or sessionStorage
   * In a real app, this would load from a file system or database
   */
  private async initializeCache(): Promise<void> {
    if (this.isInitialized) return;

    try {
      // Try to load from localStorage first (persists between sessions)
      const persistentCache = localStorage.getItem('mtb-cache');
      if (persistentCache) {
        const cacheLog: CacheLog = JSON.parse(persistentCache);
        const today = new Date().toISOString().split('T')[0];
        
        // Only load entries from today
        cacheLog.entries
          .filter(entry => entry.date === today)
          .forEach(entry => {
            const key = this.getCacheKey(entry.query.range, entry.query.subjectId);
            this.cache.set(key, entry);
          });
        
        console.log(`Loaded ${this.cache.size} cache entries from today`);
      }
    } catch (error) {
      console.warn('Failed to initialize cache from localStorage:', error);
    }

    this.isInitialized = true;
  }

  /**
   * Save cache to localStorage
   */
  private async saveCache(): Promise<void> {
    try {
      const entries = Array.from(this.cache.values());
      const cacheLog: CacheLog = { entries };
      
      localStorage.setItem('mtb-cache', JSON.stringify(cacheLog));
      
      // Also save to sessionStorage as backup
      sessionStorage.setItem('mtb-cache-session', JSON.stringify(cacheLog));
    } catch (error) {
      console.warn('Failed to save cache to localStorage:', error);
    }
  }

  /**
   * Get cached results for a query
   */
  async getCachedResults(range: string, subjectId = 'L2345'): Promise<CacheEntry | null> {
    await this.initializeCache();
    
    const key = this.getCacheKey(range, subjectId);
    const entry = this.cache.get(key);
    
    if (entry && this.isFromToday(entry)) {
      console.log(`Cache hit for ${range} (${subjectId})`);
      return entry;
    }
    
    if (entry) {
      // Remove outdated entry
      this.cache.delete(key);
      console.log(`Cache miss (outdated) for ${range} (${subjectId})`);
    } else {
      console.log(`Cache miss for ${range} (${subjectId})`);
    }
    
    return null;
  }

  /**
   * Store query results in cache
   */
  async setCachedResults(
    range: string,
    subjectId = 'L2345',
    processedVariants: ProcessedVariantCache[]
  ): Promise<void> {
    await this.initializeCache();

    const key = this.getCacheKey(range, subjectId);
    const now = new Date();

    const entry: CacheEntry = {
      timestamp: now.toISOString(),
      date: now.toISOString().split('T')[0],
      query: { range, subjectId },
      response: { processedVariants }
    };
    
    this.cache.set(key, entry);
    await this.saveCache();
    
    console.log(`Cached results for ${range} (${subjectId})`);
  }

  /**
   * Clear all cache entries
   */
  async clearCache(): Promise<void> {
    this.cache.clear();
    
    try {
      localStorage.removeItem('mtb-cache');
      sessionStorage.removeItem('mtb-cache-session');
      console.log('Cache cleared');
    } catch (error) {
      console.warn('Failed to clear cache from storage:', error);
    }
  }

  /**
   * Get cache statistics
   */
  getCacheStats(): { size: number; entries: Array<{ key: string; timestamp: string; date: string }> } {
    const entries = Array.from(this.cache.entries()).map(([key, entry]) => ({
      key,
      timestamp: entry.timestamp,
      date: entry.date
    }));
    
    return {
      size: this.cache.size,
      entries
    };
  }

  /**
   * Clean up old cache entries (older than today)
   */
  async cleanupCache(): Promise<number> {
    await this.initializeCache();
    
    const today = new Date().toISOString().split('T')[0];
    let removedCount = 0;
    
    for (const [key, entry] of this.cache.entries()) {
      if (entry.date !== today) {
        this.cache.delete(key);
        removedCount++;
      }
    }
    
    if (removedCount > 0) {
      await this.saveCache();
      console.log(`Cleaned up ${removedCount} old cache entries`);
    }
    
    return removedCount;
  }
}

export default CacheService.getInstance();