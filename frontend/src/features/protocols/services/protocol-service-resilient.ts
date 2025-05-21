/**
 * Resilient Protocol Service with circuit breaker pattern
 */

import { UUID } from 'crypto';
import { 
  ProtocolService, 
  Protocol, 
  ProtocolListParams,
  ProtocolDiff,
  ValidationResult,
  ProtocolLibraryItem,
  MockProtocolService
} from './protocol-service';
import { getApiService } from '../../api/api-service';
import { MOCK_PROTOCOLS } from '../data/mock-protocols';

/**
 * Enhanced protocol service implementation that uses the resilient API layer
 * with circuit breaker and fallback mechanisms
 */
export class ResilienceEnhancedProtocolService implements ProtocolService {
  constructor(private fallbackService?: ProtocolService) {}

  /**
   * Get appropriate API service based on circuit status
   * Uses the default implementation with backups as needed
   */
  private getApi() {
    return getApiService();
  }

  /**
   * Convert API params to URL params for the API service
   */
  private convertParams(params?: ProtocolListParams): Record<string, any> {
    if (\!params) return {};
    
    const apiParams: Record<string, any> = {};
    
    for (const [key, value] of Object.entries(params)) {
      if (value \!== undefined) {
        if (key === 'tags' && Array.isArray(value)) {
          apiParams.tags = value.join(',');
        } else {
          apiParams[key] = value;
        }
      }
    }
    
    return apiParams;
  }

  /**
   * Get protocols with pagination and filtering
   */
  async getProtocols(params?: ProtocolListParams) {
    try {
      const apiParams = this.convertParams(params);
      const response = await this.getApi().getProtocols(apiParams);
      
      return {
        data: response.data || [],
        total: response.total || 0,
        page: response.page || 1,
        per_page: response.per_page || 10,
        total_pages: response.total_pages || 1
      };
    } catch (error) {
      console.error('Error fetching protocols:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getProtocols');
        return this.fallbackService.getProtocols(params);
      }
      
      // Return empty result if no fallback
      return {
        data: [],
        total: 0,
        page: params?.page || 1,
        per_page: params?.per_page || 10,
        total_pages: 0
      };
    }
  }

  /**
   * Get protocol by ID
   */
  async getProtocol(id: UUID) {
    try {
      return await this.getApi().getProtocol(id);
    } catch (error) {
      console.error(`Error fetching protocol ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getProtocol');
        return this.fallbackService.getProtocol(id);
      }
      
      throw error;
    }
  }

  /**
   * Get versions of a protocol
   */
  async getProtocolVersions(id: UUID) {
    try {
      return await this.getApi().getProtocolVersions(id);
    } catch (error) {
      console.error(`Error fetching protocol versions for ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getProtocolVersions');
        return this.fallbackService.getProtocolVersions(id);
      }
      
      // Return empty result if no fallback
      return {
        protocol_id: id,
        versions: []
      };
    }
  }

  /**
   * Compare two protocol versions
   */
  async compareProtocolVersions(protocolId1: UUID, protocolId2: UUID) {
    try {
      return await this.getApi().compareProtocolVersions(protocolId1, protocolId2);
    } catch (error) {
      console.error(`Error comparing protocols ${protocolId1} and ${protocolId2}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for compareProtocolVersions');
        return this.fallbackService.compareProtocolVersions(protocolId1, protocolId2);
      }
      
      // Return basic empty comparison if no fallback
      return {
        id: crypto.randomUUID(),
        previous_version: 'unknown',
        current_version: 'unknown',
        changes: {
          added: [],
          removed: [],
          modified: []
        },
        compatibility_analysis: {
          compatible: true,
          breaking_changes: [],
          warnings: ['Comparison failed, using fallback data']
        },
        created_at: new Date().toISOString()
      } as ProtocolDiff;
    }
  }

  /**
   * Create new protocol
   */
  async createProtocol(protocol: Omit<Protocol, 'id'  < /dev/null |  'provenance'>) {
    try {
      return await this.getApi().createProtocol(protocol);
    } catch (error) {
      console.error('Error creating protocol:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for createProtocol');
        return this.fallbackService.createProtocol(protocol);
      }
      
      throw error;
    }
  }

  /**
   * Update existing protocol
   */
  async updateProtocol(id: UUID, protocol: Partial<Protocol>) {
    try {
      return await this.getApi().updateProtocol(id, protocol);
    } catch (error) {
      console.error(`Error updating protocol ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for updateProtocol');
        return this.fallbackService.updateProtocol(id, protocol);
      }
      
      throw error;
    }
  }

  /**
   * Delete protocol
   */
  async deleteProtocol(id: UUID) {
    try {
      await this.getApi().deleteProtocol(id);
    } catch (error) {
      console.error(`Error deleting protocol ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for deleteProtocol');
        return this.fallbackService.deleteProtocol(id);
      }
      
      throw error;
    }
  }

  /**
   * Create a new version of a protocol
   */
  async createProtocolVersion(
    protocolId: UUID,
    changes: Partial<Protocol>,
    versionNotes?: string
  ) {
    try {
      return await this.getApi().createProtocolVersion(protocolId, changes, versionNotes);
    } catch (error) {
      console.error(`Error creating protocol version for ${protocolId}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for createProtocolVersion');
        return this.fallbackService.createProtocolVersion(protocolId, changes, versionNotes);
      }
      
      throw error;
    }
  }

  /**
   * Validate a protocol
   */
  async validateProtocol(protocol: Protocol) {
    try {
      return await this.getApi().validateProtocol(protocol);
    } catch (error) {
      console.error('Error validating protocol:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for validateProtocol');
        return this.fallbackService.validateProtocol(protocol);
      }
      
      // Return basic validation if no fallback
      return {
        valid: true,
        errors: [],
        warnings: [{
          message: 'Validation service unavailable, using fallback',
          severity: 'warning' as const
        }],
        suggestions: []
      } as ValidationResult;
    }
  }

  /**
   * Export protocol in various formats
   */
  async exportProtocol(
    id: UUID,
    format: 'json' | 'yaml' | 'pdf' | 'human-readable'
  ) {
    try {
      return await this.getApi().exportProtocol(id, format);
    } catch (error) {
      console.error(`Error exporting protocol ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for exportProtocol');
        return this.fallbackService.exportProtocol(id, format);
      }
      
      // Return empty blob if no fallback
      return new Blob([`Error exporting protocol ${id}`], { type: 'text/plain' });
    }
  }

  /**
   * Import protocol from file
   */
  async importProtocol(file: File) {
    try {
      return await this.getApi().importProtocol(file);
    } catch (error) {
      console.error('Error importing protocol:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for importProtocol');
        return this.fallbackService.importProtocol(file);
      }
      
      throw error;
    }
  }

  /**
   * Get protocol templates
   */
  async getTemplates(params?: ProtocolListParams) {
    try {
      return await this.getApi().getTemplates(params);
    } catch (error) {
      console.error('Error fetching protocol templates:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getTemplates');
        return this.fallbackService.getTemplates(params);
      }
      
      // Return empty result if no fallback
      return {
        data: [],
        total: 0,
        page: params?.page || 1,
        per_page: params?.per_page || 10,
        total_pages: 0
      };
    }
  }

  /**
   * Create protocol from template
   */
  async createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>) {
    try {
      return await this.getApi().createFromTemplate(templateId, customizations);
    } catch (error) {
      console.error(`Error creating protocol from template ${templateId}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for createFromTemplate');
        return this.fallbackService.createFromTemplate(templateId, customizations);
      }
      
      throw error;
    }
  }

  /**
   * Get protocols from library
   */
  async getLibraryProtocols(params?: ProtocolListParams) {
    try {
      return await this.getApi().getLibraryProtocols(params);
    } catch (error) {
      console.error('Error fetching library protocols:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getLibraryProtocols');
        return this.fallbackService.getLibraryProtocols(params);
      }
      
      // Return empty result if no fallback
      return {
        data: [],
        total: 0,
        page: params?.page || 1,
        per_page: params?.per_page || 10,
        total_pages: 0
      };
    }
  }

  /**
   * Get protocol from library
   */
  async getLibraryProtocol(id: UUID) {
    try {
      return await this.getApi().getLibraryProtocol(id);
    } catch (error) {
      console.error(`Error fetching library protocol ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getLibraryProtocol');
        return this.fallbackService.getLibraryProtocol(id);
      }
      
      throw error;
    }
  }

  /**
   * Publish protocol to library
   */
  async publishToLibrary(
    protocolId: UUID,
    metadata: {
      category: string;
      tags: string[];
      description?: string;
    }
  ) {
    try {
      return await this.getApi().publishToLibrary(protocolId, metadata);
    } catch (error) {
      console.error(`Error publishing protocol ${protocolId} to library:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for publishToLibrary');
        return this.fallbackService.publishToLibrary(protocolId, metadata);
      }
      
      throw error;
    }
  }

  /**
   * Search protocols
   */
  async searchProtocols(query: string, params?: ProtocolListParams) {
    try {
      return await this.getApi().searchProtocols(query, params);
    } catch (error) {
      console.error('Error searching protocols:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for searchProtocols');
        return this.fallbackService.searchProtocols(query, params);
      }
      
      // Return empty result if no fallback
      return {
        data: [],
        total: 0
      };
    }
  }
}

/**
 * Create a mock service with initial data for fallback
 */
const createMockService = () => {
  const mockService = new MockProtocolService();
  
  // Pre-populate with mock data
  // @ts-ignore - Handle type incompatibility for demo
  mockService['protocols'] = MOCK_PROTOCOLS;
  
  return mockService;
};

/**
 * Create a resilient protocol service with automatic fallback
 */
export const createResilienceEnhancedProtocolService = (fallbackService?: ProtocolService) => {
  return new ResilienceEnhancedProtocolService(fallbackService || createMockService());
};
