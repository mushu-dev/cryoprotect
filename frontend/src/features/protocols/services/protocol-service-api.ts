/**
 * Protocol API Service
 * 
 * This service provides API client implementations for the Protocol service interface,
 * connecting to the backend API endpoints to manage experiment protocols.
 */

import { UUID } from 'crypto';
import { 
  Protocol, 
  ProtocolService,
  ProtocolDiff,
  ValidationResult,
  ProtocolLibraryItem,
  ProtocolListParams
} from './protocol-service';

/**
 * Implementation of the Protocol service interface using direct API calls.
 * This implementation communicates with the backend API endpoints to perform
 * protocol management operations.
 */
export class ApiProtocolService implements ProtocolService {
  private apiUrl: string;
  private headers: HeadersInit;
  
  /**
   * Creates a new instance of the API Protocol Service
   * 
   * @param baseUrl - Base URL of the API (defaults to "/api/v1")
   * @param authToken - Optional authorization token
   */
  constructor(baseUrl: string = '/api/v1', authToken?: string) {
    this.apiUrl = `${baseUrl}/protocols`;
    
    this.headers = {
      'Content-Type': 'application/json',
    };
    
    if (authToken) {
      this.headers['Authorization'] = `Bearer ${authToken}`;
    }
  }

  /**
   * Handle API responses and errors in a consistent way
   */
  private async handleResponse<T>(response: Response): Promise<T> {
    if (!response.ok) {
      const errorText = await response.text();
      let errorMessage: string;
      
      try {
        const errorData = JSON.parse(errorText);
        errorMessage = errorData.message || errorData.error || `API Error: ${response.status} ${response.statusText}`;
      } catch (e) {
        errorMessage = errorText || `API Error: ${response.status} ${response.statusText}`;
      }
      
      throw new Error(errorMessage);
    }
    
    return response.json() as Promise<T>;
  }
  
  /**
   * Fetches a list of protocols with optional filtering and pagination
   */
  async getProtocols(params?: ProtocolListParams) {
    const queryParams = new URLSearchParams();
    
    if (params) {
      Object.entries(params).forEach(([key, value]) => {
        if (value !== undefined) {
          if (Array.isArray(value)) {
            value.forEach(v => queryParams.append(key, v.toString()));
          } else {
            queryParams.append(key, value.toString());
          }
        }
      });
    }
    
    const response = await fetch(`${this.apiUrl}?${queryParams.toString()}`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      data: Protocol[];
      total: number;
      page: number;
      per_page: number;
      total_pages: number;
    }>(response);
  }
  
  /**
   * Retrieves a single protocol by ID
   */
  async getProtocol(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      headers: this.headers
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Gets all versions of a protocol
   */
  async getProtocolVersions(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}/versions`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      protocol_id: UUID;
      versions: {
        version: string;
        id: UUID;
        created_at: string;
        created_by: string;
        description?: string;
      }[];
    }>(response);
  }
  
  /**
   * Compares two protocol versions
   */
  async compareProtocolVersions(protocolId1: UUID, protocolId2: UUID) {
    const response = await fetch(`${this.apiUrl}/compare?id1=${protocolId1}&id2=${protocolId2}`, {
      headers: this.headers
    });
    
    return this.handleResponse<ProtocolDiff>(response);
  }
  
  /**
   * Creates a new protocol
   */
  async createProtocol(protocol: Omit<Protocol, 'id' | 'provenance'>) {
    const response = await fetch(this.apiUrl, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify(protocol),
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Updates an existing protocol
   */
  async updateProtocol(id: UUID, protocol: Partial<Protocol>) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      method: 'PATCH',
      headers: this.headers,
      body: JSON.stringify(protocol),
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Deletes a protocol
   */
  async deleteProtocol(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      method: 'DELETE',
      headers: this.headers
    });
    
    if (!response.ok) {
      throw new Error(`Failed to delete protocol: ${response.statusText}`);
    }
  }
  
  /**
   * Creates a new version of an existing protocol
   */
  async createProtocolVersion(protocolId: UUID, changes: Partial<Protocol>, versionNotes?: string) {
    const response = await fetch(`${this.apiUrl}/${protocolId}/versions`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify({
        changes,
        version_notes: versionNotes,
      }),
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Validates a protocol against the backend rules
   */
  async validateProtocol(protocol: Protocol) {
    const response = await fetch(`${this.apiUrl}/validate`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify(protocol),
    });
    
    return this.handleResponse<ValidationResult>(response);
  }
  
  /**
   * Exports a protocol in the requested format
   */
  async exportProtocol(id: UUID, format: 'json' | 'yaml' | 'pdf' | 'human-readable') {
    const response = await fetch(`${this.apiUrl}/${id}/export?format=${format}`, {
      headers: this.headers
    });
    
    if (!response.ok) {
      throw new Error(`Failed to export protocol: ${response.statusText}`);
    }
    
    return response.blob();
  }
  
  /**
   * Imports a protocol from a file
   */
  async importProtocol(file: File) {
    const formData = new FormData();
    formData.append('file', file);
    
    // Remove the content-type header since it will be set automatically for FormData
    const headers = { ...this.headers };
    delete headers['Content-Type'];
    
    const response = await fetch(`${this.apiUrl}/import`, {
      method: 'POST',
      headers,
      body: formData,
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Gets protocol templates
   */
  async getTemplates(params?: ProtocolListParams) {
    const queryParams = new URLSearchParams();
    queryParams.append('is_template', 'true');
    
    if (params) {
      Object.entries(params).forEach(([key, value]) => {
        if (value !== undefined && key !== 'is_template') {
          if (Array.isArray(value)) {
            value.forEach(v => queryParams.append(key, v.toString()));
          } else {
            queryParams.append(key, value.toString());
          }
        }
      });
    }
    
    const response = await fetch(`${this.apiUrl}?${queryParams.toString()}`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      data: Protocol[];
      total: number;
      page: number;
      per_page: number;
      total_pages: number;
    }>(response);
  }
  
  /**
   * Creates a new protocol from a template
   */
  async createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>) {
    const response = await fetch(`${this.apiUrl}/create-from-template`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify({
        template_id: templateId,
        customizations,
      }),
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Gets protocols from the library
   */
  async getLibraryProtocols(params?: ProtocolListParams) {
    const queryParams = new URLSearchParams();
    
    if (params) {
      Object.entries(params).forEach(([key, value]) => {
        if (value !== undefined) {
          if (Array.isArray(value)) {
            value.forEach(v => queryParams.append(key, v.toString()));
          } else {
            queryParams.append(key, value.toString());
          }
        }
      });
    }
    
    const response = await fetch(`${this.apiUrl}/library?${queryParams.toString()}`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      data: ProtocolLibraryItem[];
      total: number;
      page: number;
      per_page: number;
      total_pages: number;
    }>(response);
  }
  
  /**
   * Gets a protocol from the library
   */
  async getLibraryProtocol(id: UUID) {
    const response = await fetch(`${this.apiUrl}/library/${id}`, {
      headers: this.headers
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Publishes a protocol to the library
   */
  async publishToLibrary(protocolId: UUID, metadata: {
    category: string;
    tags: string[];
    description?: string;
  }) {
    const response = await fetch(`${this.apiUrl}/${protocolId}/publish`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify(metadata),
    });
    
    return this.handleResponse<ProtocolLibraryItem>(response);
  }
  
  /**
   * Searches protocols
   */
  async searchProtocols(query: string, params?: ProtocolListParams) {
    const queryParams = new URLSearchParams();
    queryParams.append('query', query);
    
    if (params) {
      Object.entries(params).forEach(([key, value]) => {
        if (value !== undefined) {
          if (Array.isArray(value)) {
            value.forEach(v => queryParams.append(key, v.toString()));
          } else {
            queryParams.append(key, value.toString());
          }
        }
      });
    }
    
    const response = await fetch(`${this.apiUrl}/search?${queryParams.toString()}`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      data: Protocol[];
      total: number;
    }>(response);
  }
  
  /**
   * Gets protocols for a specific mixture
   */
  async getProtocolsForMixture(mixtureId: UUID) {
    const response = await fetch(`${this.apiUrl}/mixtures/${mixtureId}/protocols`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      mixture_id: UUID;
      mixture_name: string;
      protocols: Protocol[];
    }>(response);
  }
  
  /**
   * Gets sensitivity profiles for different sample types
   */
  async getSensitivityProfiles() {
    const response = await fetch(`${this.apiUrl}/sensitivity-profiles`, {
      headers: this.headers
    });
    
    return this.handleResponse<{
      profiles: Record<string, {
        name: string;
        description: string;
        osmotic_tolerance: number;
        max_step_size: number;
        time_per_step: number;
        cooling_rate: number;
        warming_rate: number;
        notes: string;
      }>;
    }>(response);
  }
  
  /**
   * Designs a concentration gradient protocol for a mixture
   */
  async designProtocol(
    mixtureId: UUID, 
    parameters: {
      target_concentration: number;
      sample_type: string;
      starting_temperature: number;
      target_temperature?: number;
      step_count?: number;
      custom_sensitivity?: Record<string, any>;
    }
  ) {
    const response = await fetch(`${this.apiUrl}/design/${mixtureId}`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify(parameters),
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Saves a protocol for a mixture
   */
  async saveProtocol(
    mixtureId: UUID,
    protocol: {
      name: string;
      description?: string;
      target_concentration: number;
      sample_type: string;
      starting_temperature: number;
      target_temperature?: number;
      step_count?: number;
      steps?: any[];
      custom_sensitivity?: Record<string, any>;
    }
  ) {
    const response = await fetch(`${this.apiUrl}/save/${mixtureId}`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify(protocol),
    });
    
    return this.handleResponse<Protocol>(response);
  }
  
  /**
   * Compares multiple protocols
   */
  async compareProtocols(protocolIds: UUID[]) {
    const response = await fetch(`${this.apiUrl}/compare`, {
      method: 'POST',
      headers: this.headers,
      body: JSON.stringify({ protocol_ids: protocolIds }),
    });
    
    return this.handleResponse<{
      protocols: Protocol[];
      comparison: Record<string, any>;
    }>(response);
  }
}