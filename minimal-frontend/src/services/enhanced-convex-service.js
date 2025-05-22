/**
 * Enhanced Convex Service for CryoProtect
 * 
 * Provides direct connection to the enhanced Convex database with:
 * - Real-time molecular properties and cryoprotectant scores
 * - Advanced caching integration
 * - Quality validation and monitoring
 * - Performance optimizations
 */

// Base URL for Convex API
const CONVEX_URL = process.env.NEXT_PUBLIC_CONVEX_URL || 'https://upbeat-parrot-866.convex.cloud';

/**
 * Enhanced Convex Client for direct API communication
 */
class EnhancedConvexClient {
  constructor() {
    this.baseUrl = CONVEX_URL;
    this.cache = new Map();
    this.cacheTimeout = 5 * 60 * 1000; // 5 minutes
    
    // Performance metrics
    this.metrics = {
      requests: 0,
      cacheHits: 0,
      cacheMisses: 0,
      avgResponseTime: 0
    };
  }

  /**
   * Call a Convex function via HTTP API
   */
  async callFunction(functionPath, args = {}, isQuery = true) {
    const startTime = Date.now();
    this.metrics.requests++;

    // Check cache first for queries
    if (isQuery) {
      const cacheKey = this._generateCacheKey(functionPath, args);
      const cached = this._getFromCache(cacheKey);
      if (cached) {
        this.metrics.cacheHits++;
        return cached;
      }
      this.metrics.cacheMisses++;
    }

    try {
      const endpoint = isQuery ? 'query' : 'mutation';
      const response = await fetch(`${this.baseUrl}/api/${endpoint}`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          path: functionPath,
          args: args,
          format: 'json'
        }),
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const result = await response.json();
      
      // Handle Convex response format
      const data = result.status === 'success' ? result.value : result;
      
      // Cache successful query results
      if (isQuery && data) {
        const cacheKey = this._generateCacheKey(functionPath, args);
        this._setToCache(cacheKey, data);
      }

      // Update performance metrics
      const responseTime = Date.now() - startTime;
      this._updateResponseTime(responseTime);

      return data;
    } catch (error) {
      console.error(`Convex ${isQuery ? 'query' : 'mutation'} failed:`, error);
      throw error;
    }
  }

  /**
   * Generate cache key for function call
   */
  _generateCacheKey(functionPath, args) {
    return `${functionPath}:${JSON.stringify(args)}`;
  }

  /**
   * Get value from cache if not expired
   */
  _getFromCache(key) {
    const entry = this.cache.get(key);
    if (entry && Date.now() - entry.timestamp < this.cacheTimeout) {
      return entry.data;
    }
    if (entry) {
      this.cache.delete(key);
    }
    return null;
  }

  /**
   * Set value to cache with timestamp
   */
  _setToCache(key, data) {
    this.cache.set(key, {
      data,
      timestamp: Date.now()
    });
  }

  /**
   * Update average response time metric
   */
  _updateResponseTime(responseTime) {
    this.metrics.avgResponseTime = (
      (this.metrics.avgResponseTime * (this.metrics.requests - 1) + responseTime) / 
      this.metrics.requests
    );
  }

  /**
   * Clear cache
   */
  clearCache() {
    this.cache.clear();
  }

  /**
   * Get performance metrics
   */
  getMetrics() {
    const hitRate = this.metrics.requests > 0 
      ? (this.metrics.cacheHits / this.metrics.requests * 100).toFixed(1)
      : 0;
    
    return {
      ...this.metrics,
      hitRate: parseFloat(hitRate),
      cacheSize: this.cache.size
    };
  }
}

// Create global client instance
const convexClient = new EnhancedConvexClient();

/**
 * Enhanced Molecule Service with Convex Integration
 */
export const enhancedConvexMoleculeService = {
  /**
   * Get all molecules with enhanced properties and scores
   */
  async getMolecules(params = {}) {
    try {
      const molecules = await convexClient.callFunction('molecules:listMolecules', {
        limit: params.limit || 20
      });

      if (!molecules) {
        return {
          molecules: [],
          total: 0,
          page: params.page || 1,
          limit: params.limit || 20
        };
      }

      // Get enhanced data for each molecule in parallel
      const enhancedMolecules = await Promise.all(
        molecules.map(async (molecule) => {
          try {
            const enhanced = await convexClient.callFunction('molecules:getMoleculeWithProperties', {
              id: molecule._id
            });

            return this._formatMoleculeForFrontend(enhanced || molecule);
          } catch (error) {
            console.warn(`Failed to get enhanced data for molecule ${molecule._id}:`, error);
            return this._formatMoleculeForFrontend(molecule);
          }
        })
      );

      // Apply search filter client-side if provided
      let filteredMolecules = enhancedMolecules;
      if (params.search) {
        const searchTerm = params.search.toLowerCase();
        filteredMolecules = enhancedMolecules.filter(molecule =>
          molecule.name.toLowerCase().includes(searchTerm) ||
          (molecule.formula && molecule.formula.toLowerCase().includes(searchTerm)) ||
          (molecule.smiles && molecule.smiles.toLowerCase().includes(searchTerm))
        );
      }

      // Apply cryoprotectant filter if requested
      if (params.is_cryoprotectant) {
        filteredMolecules = filteredMolecules.filter(molecule =>
          molecule.cryoprotectantScore && molecule.cryoprotectantScore > 0
        );
      }

      // Sort molecules
      if (params.sort_by) {
        this._sortMolecules(filteredMolecules, params.sort_by, params.sort_order || 'asc');
      }

      // Apply pagination
      const page = params.page || 1;
      const limit = params.limit || 20;
      const startIndex = (page - 1) * limit;
      const paginatedMolecules = filteredMolecules.slice(startIndex, startIndex + limit);

      return {
        molecules: paginatedMolecules,
        total: filteredMolecules.length,
        page: page,
        limit: limit
      };
    } catch (error) {
      console.error('Error fetching enhanced molecules:', error);
      throw error;
    }
  },

  /**
   * Get a single molecule with all enhanced data
   */
  async getMolecule(id) {
    try {
      const molecule = await convexClient.callFunction('molecules:getMoleculeWithProperties', {
        id: id
      });

      if (!molecule) {
        throw new Error('Molecule not found');
      }

      return this._formatMoleculeForFrontend(molecule);
    } catch (error) {
      console.error(`Error fetching enhanced molecule ${id}:`, error);
      throw error;
    }
  },

  /**
   * Get molecule properties
   */
  async getMoleculeProperties(id) {
    try {
      const molecule = await this.getMolecule(id);
      return {
        properties: molecule.properties || []
      };
    } catch (error) {
      console.error(`Error fetching properties for molecule ${id}:`, error);
      throw error;
    }
  },

  /**
   * Search molecules with enhanced scoring
   */
  async searchMolecules(query, limit = 10) {
    try {
      const allMolecules = await this.getMolecules({ limit: 1000 });
      
      const searchTerm = query.toLowerCase();
      const results = allMolecules.molecules.filter(molecule =>
        molecule.name.toLowerCase().includes(searchTerm) ||
        (molecule.formula && molecule.formula.toLowerCase().includes(searchTerm)) ||
        (molecule.smiles && molecule.smiles.toLowerCase().includes(searchTerm))
      );

      // Sort by relevance (name matches first, then by cryoprotectant score)
      results.sort((a, b) => {
        const aNameMatch = a.name.toLowerCase().includes(searchTerm) ? 1 : 0;
        const bNameMatch = b.name.toLowerCase().includes(searchTerm) ? 1 : 0;
        
        if (aNameMatch !== bNameMatch) {
          return bNameMatch - aNameMatch;
        }
        
        return (b.cryoprotectantScore || 0) - (a.cryoprotectantScore || 0);
      });

      return {
        results: results.slice(0, limit)
      };
    } catch (error) {
      console.error('Error searching enhanced molecules:', error);
      throw error;
    }
  },

  /**
   * Get molecules ranked by cryoprotectant effectiveness
   */
  async getCryoprotectantRanking(category = null, limit = 10) {
    try {
      const molecules = await convexClient.callFunction('molecules:listMoleculesWithScores', {
        category: category,
        limit: limit
      });

      return {
        molecules: molecules.map(mol => this._formatMoleculeForFrontend(mol))
      };
    } catch (error) {
      console.error('Error fetching cryoprotectant ranking:', error);
      throw error;
    }
  },

  /**
   * Get quality metrics for the database
   */
  async getQualityMetrics() {
    try {
      // This would call a quality assessment function if available
      const molecules = await convexClient.callFunction('molecules:listMolecules', {
        limit: 1000
      });

      let totalMolecules = 0;
      let moleculesWithProperties = 0;
      let moleculesWithScores = 0;
      let avgQualityScore = 0;

      if (molecules) {
        totalMolecules = molecules.length;
        
        for (const molecule of molecules) {
          const enhanced = await convexClient.callFunction('molecules:getMoleculeWithProperties', {
            id: molecule._id
          });

          if (enhanced?.properties) {
            moleculesWithProperties++;
          }
          if (enhanced?.scores) {
            moleculesWithScores++;
            avgQualityScore += enhanced.scores.overallScore || 0;
          }
        }
        
        avgQualityScore = moleculesWithScores > 0 ? avgQualityScore / moleculesWithScores : 0;
      }

      return {
        totalMolecules,
        moleculesWithProperties,
        moleculesWithScores,
        avgQualityScore: Math.round(avgQualityScore * 100) / 100,
        propertiesCompletionRate: totalMolecules > 0 ? (moleculesWithProperties / totalMolecules * 100) : 0,
        scoresCompletionRate: totalMolecules > 0 ? (moleculesWithScores / totalMolecules * 100) : 0
      };
    } catch (error) {
      console.error('Error fetching quality metrics:', error);
      throw error;
    }
  },

  /**
   * Format molecule data for frontend consumption
   */
  _formatMoleculeForFrontend(molecule) {
    if (!molecule) return null;

    // Extract properties
    const properties = molecule.properties || {};
    const scores = molecule.scores || {};

    return {
      id: molecule._id,
      name: molecule.name,
      formula: molecule.formula,
      smiles: molecule.canonicalSmiles,
      inchiKey: molecule.inchiKey,
      pubchemCid: molecule.pubchemCid,
      status: molecule.status,
      createdAt: molecule.createdAt,
      updatedAt: molecule.updatedAt,
      
      // Molecular properties
      molecularWeight: properties.molecularWeight,
      exactMass: properties.exactMass,
      logP: properties.logP,
      tpsa: properties.tpsa,
      hbondDonors: properties.hbondDonors,
      hbondAcceptors: properties.hbondAcceptors,
      rotatableBonds: properties.rotatableBonds,
      aromaticRings: properties.aromaticRings,
      aliphaticRings: properties.aliphaticRings,
      complexity: properties.complexity,
      heavyAtomCount: properties.heavyAtomCount,
      morganFingerprint: properties.morganFingerprint,
      rdkitFingerprint: properties.rdkitFingerprint,
      
      // Cryoprotectant scores
      cryoprotectantScore: scores.overallScore,
      cryoprotectantCategory: scores.category,
      glassTempScore: scores.glassTempScore,
      viscosityScore: scores.viscosityScore,
      permeabilityScore: scores.permeabilityScore,
      toxicityScore: scores.toxicityScore,
      confidence: scores.confidence,
      
      // Combined properties for easy access
      properties: {
        ...properties,
        calculatedAt: properties.calculatedAt,
        calculationVersion: properties.calculationVersion
      },
      
      cryoprotectantScores: {
        ...scores,
        calculatedAt: scores.calculatedAt,
        scoringAlgorithmVersion: scores.scoringAlgorithmVersion
      }
    };
  },

  /**
   * Sort molecules by specified field
   */
  _sortMolecules(molecules, sortBy, sortOrder) {
    const direction = sortOrder === 'desc' ? -1 : 1;
    
    molecules.sort((a, b) => {
      let aValue = a[sortBy];
      let bValue = b[sortBy];
      
      // Handle special cases
      if (sortBy === 'molecular_weight') {
        aValue = a.molecularWeight;
        bValue = b.molecularWeight;
      } else if (sortBy === 'created_at') {
        aValue = new Date(a.createdAt);
        bValue = new Date(b.createdAt);
      }
      
      // Handle null/undefined values
      if (aValue == null && bValue == null) return 0;
      if (aValue == null) return 1;
      if (bValue == null) return -1;
      
      // Compare values
      if (typeof aValue === 'string') {
        return direction * aValue.localeCompare(bValue);
      } else if (typeof aValue === 'number' || aValue instanceof Date) {
        return direction * (aValue - bValue);
      }
      
      return 0;
    });
  },

  /**
   * Get client performance metrics
   */
  getClientMetrics() {
    return convexClient.getMetrics();
  },

  /**
   * Clear client cache
   */
  clearCache() {
    convexClient.clearCache();
  }
};

/**
 * Initialize enhanced Convex service
 */
export const initializeEnhancedConvexService = () => {
  console.log('🚀 Enhanced Convex service initialized');
  console.log(`📡 Connected to: ${CONVEX_URL}`);
  
  // Test connection
  enhancedConvexMoleculeService.getMolecules({ limit: 1 })
    .then(() => {
      console.log('✅ Enhanced Convex connection verified');
    })
    .catch((error) => {
      console.warn('⚠️ Enhanced Convex connection failed:', error);
    });
};

// Auto-initialize on import
if (typeof window !== 'undefined') {
  initializeEnhancedConvexService();
}