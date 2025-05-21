/**
 * Protocol Service with Context Integration
 * 
 * This service implements the Protocol service interface using the application's
 * resilient API client from the API context.
 */

import { UUID } from 'crypto';
import { ResilientApiClient } from '@/services/resilient-api-client';
import { 
  Protocol, 
  ProtocolService,
  ProtocolDiff,
  ValidationResult,
  ProtocolLibraryItem,
  ProtocolListParams
} from './protocol-service';

/**
 * Implementation of the Protocol service using the context-provided resilient API client
 */
export class ContextProtocolService implements ProtocolService {
  private apiClient: ResilientApiClient;
  private apiBaseUrl: string;
  
  /**
   * Creates a new instance of the Context Protocol Service
   * 
   * @param apiClient - The resilient API client from context
   * @param basePathPrefix - Optional base path prefix (defaults to "/protocols")
   */
  constructor(apiClient: ResilientApiClient, basePathPrefix: string = '/protocols') {
    this.apiClient = apiClient;
    this.apiBaseUrl = basePathPrefix;
  }
  
  /**
   * Helper method to build query params
   */
  private buildQueryParams(params?: ProtocolListParams): string {
    if (!params) return '';
    
    const queryParams = new URLSearchParams();
    
    Object.entries(params).forEach(([key, value]) => {
      if (value !== undefined) {
        if (Array.isArray(value)) {
          value.forEach(v => queryParams.append(key, v.toString()));
        } else {
          queryParams.append(key, value.toString());
        }
      }
    });
    
    const queryString = queryParams.toString();
    return queryString ? `?${queryString}` : '';
  }
  
  /**
   * Fetches a list of protocols with optional filtering and pagination
   */
  async getProtocols(params?: ProtocolListParams) {
    const queryString = this.buildQueryParams(params);
    return this.apiClient.get<{
      data: Protocol[];
      total: number;
      page: number;
      per_page: number;
      total_pages: number;
    }>(`${this.apiBaseUrl}${queryString}`);
  }
  
  /**
   * Retrieves a single protocol by ID
   */
  async getProtocol(id: UUID) {
    return this.apiClient.get<Protocol>(`${this.apiBaseUrl}/${id}`);
  }
  
  /**
   * Gets all versions of a protocol
   */
  async getProtocolVersions(id: UUID) {
    return this.apiClient.get<{
      protocol_id: UUID;
      versions: {
        version: string;
        id: UUID;
        created_at: string;
        created_by: string;
        description?: string;
      }[];
    }>(`${this.apiBaseUrl}/${id}/versions`);
  }
  
  /**
   * Compares two protocol versions
   */
  async compareProtocolVersions(protocolId1: UUID, protocolId2: UUID) {
    return this.apiClient.get<ProtocolDiff>(
      `${this.apiBaseUrl}/compare?id1=${protocolId1}&id2=${protocolId2}`
    );
  }
  
  /**
   * Creates a new protocol
   */
  async createProtocol(protocol: Omit<Protocol, 'id' | 'provenance'>) {
    return this.apiClient.post<Protocol>(this.apiBaseUrl, protocol);
  }
  
  /**
   * Updates an existing protocol
   */
  async updateProtocol(id: UUID, protocol: Partial<Protocol>) {
    return this.apiClient.patch<Protocol>(`${this.apiBaseUrl}/${id}`, protocol);
  }
  
  /**
   * Deletes a protocol
   */
  async deleteProtocol(id: UUID) {
    await this.apiClient.delete(`${this.apiBaseUrl}/${id}`);
  }
  
  /**
   * Creates a new version of an existing protocol
   */
  async createProtocolVersion(protocolId: UUID, changes: Partial<Protocol>, versionNotes?: string) {
    return this.apiClient.post<Protocol>(`${this.apiBaseUrl}/${protocolId}/versions`, {
      changes,
      version_notes: versionNotes,
    });
  }
  
  /**
   * Validates a protocol against the backend rules
   */
  async validateProtocol(protocol: Protocol) {
    return this.apiClient.post<ValidationResult>(`${this.apiBaseUrl}/validate`, protocol);
  }
  
  /**
   * Exports a protocol in the requested format
   */
  async exportProtocol(id: UUID, format: 'json' | 'yaml' | 'pdf' | 'human-readable') {
    const response = await fetch(`${this.apiBaseUrl}/${id}/export?format=${format}`, {
      headers: {
        Authorization: `Bearer ${localStorage.getItem('auth_token') || ''}`,
      },
    });
    
    if (!response.ok) {
      throw new Error(`Failed to export protocol: ${response.statusText}`);
    }
    
    return await response.blob();
  }
  
  /**
   * Imports a protocol from a file
   */
  async importProtocol(file: File) {
    const formData = new FormData();
    formData.append('file', file);
    
    const response = await fetch(`${this.apiBaseUrl}/import`, {
      method: 'POST',
      headers: {
        Authorization: `Bearer ${localStorage.getItem('auth_token') || ''}`,
      },
      body: formData,
    });
    
    if (!response.ok) {
      throw new Error(`Failed to import protocol: ${response.statusText}`);
    }
    
    return response.json();
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
    
    return this.apiClient.get<{
      data: Protocol[];
      total: number;
      page: number;
      per_page: number;
      total_pages: number;
    }>(`${this.apiBaseUrl}?${queryParams.toString()}`);
  }
  
  /**
   * Creates a new protocol from a template
   */
  async createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>) {
    return this.apiClient.post<Protocol>(`${this.apiBaseUrl}/create-from-template`, {
      template_id: templateId,
      customizations,
    });
  }
  
  /**
   * Gets protocols from the library
   */
  async getLibraryProtocols(params?: ProtocolListParams) {
    const queryString = this.buildQueryParams(params);
    return this.apiClient.get<{
      data: ProtocolLibraryItem[];
      total: number;
      page: number;
      per_page: number;
      total_pages: number;
    }>(`${this.apiBaseUrl}/library${queryString}`);
  }
  
  /**
   * Gets a protocol from the library
   */
  async getLibraryProtocol(id: UUID) {
    return this.apiClient.get<Protocol>(`${this.apiBaseUrl}/library/${id}`);
  }
  
  /**
   * Publishes a protocol to the library
   */
  async publishToLibrary(protocolId: UUID, metadata: {
    category: string;
    tags: string[];
    description?: string;
  }) {
    return this.apiClient.post<ProtocolLibraryItem>(
      `${this.apiBaseUrl}/${protocolId}/publish`, 
      metadata
    );
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
    
    return this.apiClient.get<{
      data: Protocol[];
      total: number;
    }>(`${this.apiBaseUrl}/search?${queryParams.toString()}`);
  }
}