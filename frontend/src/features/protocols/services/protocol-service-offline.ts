/**
 * Offline-Capable Protocol Service
 * 
 * This service implements the Protocol service interface using the sync service
 * for offline capability, allowing the application to work even when offline.
 */

import { UUID } from 'crypto';
import { v4 as uuidv4 } from 'uuid';
import { SyncService } from '@/services/sync-service';
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
 * Implementation of the Protocol service using the sync service for offline capability
 */
export class OfflineProtocolService implements ProtocolService {
  private apiClient: ResilientApiClient;
  private syncService: SyncService;
  private entityType = 'protocols';
  
  /**
   * Creates a new instance of the Offline Protocol Service
   * 
   * @param apiClient - The resilient API client
   * @param syncService - The sync service for offline capability
   */
  constructor(apiClient: ResilientApiClient, syncService: SyncService) {
    this.apiClient = apiClient;
    this.syncService = syncService;
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
    try {
      // Try to get from API first
      const queryString = this.buildQueryParams(params);
      const response = await this.apiClient.get<{
        data: Protocol[];
        total: number;
        page: number;
        per_page: number;
        total_pages: number;
      }>(`/${this.entityType}${queryString}`);
      
      // Cache each protocol individually
      for (const protocol of response.data) {
        await this.syncService.saveEntity(this.entityType, protocol.id, protocol);
      }
      
      return response;
    } catch (error) {
      // If API fails, try to get from offline storage
      console.log('Using offline data for protocols list');
      
      const protocols = await this.syncService.getAllEntities<Protocol>(this.entityType);
      
      // Apply filters if necessary
      let filteredProtocols = [...protocols];
      
      if (params) {
        // Filter by experiment type
        if (params.experiment_type) {
          filteredProtocols = filteredProtocols.filter(p => 
            p.steps.some(s => s.parameters?.experiment_type === params.experiment_type)
          );
        }
        
        // Filter by template status
        if (params.is_template !== undefined) {
          filteredProtocols = filteredProtocols.filter(p => 
            p.is_template === params.is_template
          );
        }
        
        // Filter by author
        if (params.author) {
          filteredProtocols = filteredProtocols.filter(p => 
            p.created_by === params.author
          );
        }
        
        // Filter by tags
        if (params.tags && params.tags.length > 0) {
          filteredProtocols = filteredProtocols.filter(p => 
            p.tags && params.tags?.some(tag => p.tags?.includes(tag))
          );
        }
        
        // Filter by search
        if (params.search) {
          const search = params.search.toLowerCase();
          filteredProtocols = filteredProtocols.filter(p => 
            p.name.toLowerCase().includes(search) || 
            (p.description && p.description.toLowerCase().includes(search)) ||
            p.steps.some(s => 
              s.name.toLowerCase().includes(search) ||
              (s.description && s.description.toLowerCase().includes(search))
            )
          );
        }
        
        // Sort if needed
        if (params.sort_by) {
          const sortBy = params.sort_by as keyof Protocol;
          const sortOrder = params.sort_order === 'desc' ? -1 : 1;
          
          filteredProtocols.sort((a, b) => {
            const valueA = a[sortBy];
            const valueB = b[sortBy];
            
            if (typeof valueA === 'string' && typeof valueB === 'string') {
              return sortOrder * valueA.localeCompare(valueB);
            }
            
            if (valueA < valueB) return -1 * sortOrder;
            if (valueA > valueB) return 1 * sortOrder;
            return 0;
          });
        }
      }
      
      // Apply pagination
      const page = params?.page || 1;
      const perPage = params?.per_page || 10;
      const startIndex = (page - 1) * perPage;
      const endIndex = startIndex + perPage;
      
      const paginatedProtocols = filteredProtocols.slice(startIndex, endIndex);
      
      return {
        data: paginatedProtocols,
        total: filteredProtocols.length,
        page,
        per_page: perPage,
        total_pages: Math.ceil(filteredProtocols.length / perPage)
      };
    }
  }
  
  /**
   * Retrieves a single protocol by ID
   */
  async getProtocol(id: UUID) {
    try {
      // Try to get from API first
      const response = await this.apiClient.get<Protocol>(`/${this.entityType}/${id}`);
      
      // Cache the protocol
      await this.syncService.saveEntity(this.entityType, id, response);
      
      return response;
    } catch (error) {
      // If API fails, try to get from offline storage
      console.log(`Using offline data for protocol ${id}`);
      
      const protocol = await this.syncService.getEntity<Protocol>(this.entityType, id);
      
      if (!protocol) {
        throw new Error(`Protocol with ID ${id} not found`);
      }
      
      return protocol;
    }
  }
  
  /**
   * Gets all versions of a protocol
   */
  async getProtocolVersions(id: UUID) {
    try {
      // This is a read-only operation that requires server access
      return await this.apiClient.get<{
        protocol_id: UUID;
        versions: {
          version: string;
          id: UUID;
          created_at: string;
          created_by: string;
          description?: string;
        }[];
      }>(`/${this.entityType}/${id}/versions`);
    } catch (error) {
      // For now, we can't provide this feature offline
      console.error('Cannot get protocol versions while offline', error);
      
      // Return a minimal response with just the current protocol
      const protocol = await this.syncService.getEntity<Protocol>(this.entityType, id);
      
      if (!protocol) {
        throw new Error(`Protocol with ID ${id} not found`);
      }
      
      return {
        protocol_id: id,
        versions: [
          {
            version: protocol.version,
            id: protocol.id,
            created_at: protocol.created_at,
            created_by: protocol.created_by,
            description: protocol.description,
          }
        ]
      };
    }
  }
  
  /**
   * Compares two protocol versions
   */
  async compareProtocolVersions(protocolId1: UUID, protocolId2: UUID) {
    try {
      // This is a read-only operation that requires server access
      return await this.apiClient.get<ProtocolDiff>(
        `/${this.entityType}/compare?id1=${protocolId1}&id2=${protocolId2}`
      );
    } catch (error) {
      // For now, we can't provide this feature offline
      console.error('Cannot compare protocol versions while offline', error);
      
      // Return a minimal diff indicating the feature is unavailable
      return {
        id: uuidv4() as UUID,
        previous_version: 'unavailable',
        current_version: 'unavailable',
        changes: {
          added: [],
          removed: [],
          modified: [],
        },
        compatibility_analysis: {
          compatible: true,
          breaking_changes: [],
          warnings: ['Version comparison is unavailable while offline'],
        },
        created_at: new Date().toISOString(),
      };
    }
  }
  
  /**
   * Creates a new protocol
   */
  async createProtocol(protocol: Omit<Protocol, 'id' | 'provenance'>) {
    try {
      // Try to create on the server first
      const response = await this.apiClient.post<Protocol>(`/${this.entityType}`, protocol);
      
      // Cache the created protocol
      await this.syncService.saveEntity(this.entityType, response.id, response);
      
      return response;
    } catch (error) {
      // If API fails, create locally with a temporary ID
      console.log('Creating protocol offline');
      
      const tempId = uuidv4() as UUID;
      const timestamp = new Date().toISOString();
      const deviceId = this.syncService.getDeviceId();
      
      const newProtocol: Protocol = {
        ...protocol as any,
        id: tempId,
        version: protocol.version || '1.0.0',
        created_at: timestamp,
        updated_at: timestamp,
        provenance: {
          created_by: protocol.created_by,
          created_at: timestamp,
          source: 'offline-creation',
          method: 'ui',
          version: '1.0',
          device_id: deviceId,
        },
        _temp_id: true, // Mark as temporary
      };
      
      // Save to offline storage
      await this.syncService.saveEntity(this.entityType, tempId, newProtocol);
      
      return newProtocol;
    }
  }
  
  /**
   * Updates an existing protocol
   */
  async updateProtocol(id: UUID, protocol: Partial<Protocol>) {
    try {
      // Try to update on the server first
      const response = await this.apiClient.patch<Protocol>(`/${this.entityType}/${id}`, protocol);
      
      // Cache the updated protocol
      await this.syncService.saveEntity(this.entityType, id, response);
      
      return response;
    } catch (error) {
      // If API fails, update locally
      console.log(`Updating protocol ${id} offline`);
      
      // Get the current protocol
      const existingProtocol = await this.syncService.getEntity<Protocol>(this.entityType, id);
      
      if (!existingProtocol) {
        throw new Error(`Protocol with ID ${id} not found`);
      }
      
      // Update the protocol
      const updatedProtocol: Protocol = {
        ...existingProtocol,
        ...protocol,
        updated_at: new Date().toISOString(),
      };
      
      // Save to offline storage
      await this.syncService.saveEntity(this.entityType, id, updatedProtocol);
      
      return updatedProtocol;
    }
  }
  
  /**
   * Deletes a protocol
   */
  async deleteProtocol(id: UUID) {
    try {
      // Try to delete on the server first
      await this.apiClient.delete(`/${this.entityType}/${id}`);
      
      // Delete from offline storage
      await this.syncService.deleteEntity(this.entityType, id);
    } catch (error) {
      // If API fails, mark for deletion
      console.log(`Deleting protocol ${id} offline`);
      
      // Delete from offline storage and queue for server deletion when online
      await this.syncService.deleteEntity(this.entityType, id);
    }
  }
  
  /**
   * Creates a new version of an existing protocol
   */
  async createProtocolVersion(protocolId: UUID, changes: Partial<Protocol>, versionNotes?: string) {
    try {
      // Try to create on the server first
      const response = await this.apiClient.post<Protocol>(
        `/${this.entityType}/${protocolId}/versions`,
        {
          changes,
          version_notes: versionNotes,
        }
      );
      
      // Cache the new version
      await this.syncService.saveEntity(this.entityType, response.id, response);
      
      return response;
    } catch (error) {
      // If API fails, create locally with a temporary ID
      console.log(`Creating protocol version for ${protocolId} offline`);
      
      // Get the current protocol
      const existingProtocol = await this.syncService.getEntity<Protocol>(this.entityType, protocolId);
      
      if (!existingProtocol) {
        throw new Error(`Protocol with ID ${protocolId} not found`);
      }
      
      // Parse current version and increment
      const versionParts = existingProtocol.version.split('.');
      const newVersion = `${versionParts[0]}.${parseInt(versionParts[1], 10) + 1}.0`;
      
      const tempId = uuidv4() as UUID;
      const timestamp = new Date().toISOString();
      const deviceId = this.syncService.getDeviceId();
      
      // Create the new version
      const newVersionProtocol: Protocol = {
        ...existingProtocol,
        ...changes,
        id: tempId,
        version: newVersion,
        parent_version_id: protocolId,
        created_at: timestamp,
        updated_at: timestamp,
        provenance: {
          created_by: existingProtocol.created_by,
          created_at: timestamp,
          source: 'offline-version-creation',
          method: 'ui',
          version: '1.0',
          device_id: deviceId,
          version_notes: versionNotes,
        },
        _temp_id: true, // Mark as temporary
      };
      
      // Save to offline storage
      await this.syncService.saveEntity(this.entityType, tempId, newVersionProtocol);
      
      return newVersionProtocol;
    }
  }
  
  /**
   * Validates a protocol against the backend rules
   */
  async validateProtocol(protocol: Protocol) {
    try {
      // This operation requires server access for proper validation
      return await this.apiClient.post<ValidationResult>(
        `/${this.entityType}/validate`,
        protocol
      );
    } catch (error) {
      // Perform basic validation offline
      console.log('Validating protocol offline with basic checks');
      
      const errors: any[] = [];
      const warnings: any[] = [];
      const suggestions: any[] = [];
      
      // Check for required fields
      if (!protocol.name || protocol.name.trim() === '') {
        errors.push({
          field: 'name',
          message: 'Protocol name is required',
          severity: 'error',
        });
      }
      
      // Check steps
      if (!protocol.steps || protocol.steps.length === 0) {
        errors.push({
          message: 'Protocol must have at least one step',
          severity: 'error',
        });
      } else {
        // Check each step
        protocol.steps.forEach((step, index) => {
          if (!step.name || step.name.trim() === '') {
            errors.push({
              step_id: step.id,
              field: 'name',
              message: `Step ${index + 1} name is required`,
              severity: 'error',
            });
          }
          
          // Add a warning if no description
          if (!step.description || step.description.trim() === '') {
            warnings.push({
              step_id: step.id,
              field: 'description',
              message: `Step ${index + 1} should have a description`,
            });
          }
          
          // Check required parameters
          if (!step.duration) {
            warnings.push({
              step_id: step.id,
              field: 'duration',
              message: `Step ${index + 1} should specify a duration`,
            });
          }
        });
      }
      
      // Add a suggestion for longer protocols
      if (protocol.steps && protocol.steps.length < 3) {
        suggestions.push({
          message: 'Consider adding more steps for a more detailed protocol',
        });
      }
      
      return {
        valid: errors.length === 0,
        errors,
        warnings,
        suggestions,
      };
    }
  }
  
  /**
   * Exports a protocol in the requested format
   */
  async exportProtocol(id: UUID, format: 'json' | 'yaml' | 'pdf' | 'human-readable') {
    try {
      // This operation requires server access for proper formatting
      const response = await fetch(`${this.apiClient.baseURL}/${this.entityType}/${id}/export?format=${format}`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('auth_token') || ''}`,
        },
      });
      
      if (!response.ok) {
        throw new Error(`Failed to export protocol: ${response.statusText}`);
      }
      
      return await response.blob();
    } catch (error) {
      // For JSON format, we can provide an offline alternative
      if (format === 'json') {
        console.log(`Exporting protocol ${id} as JSON offline`);
        
        // Get the protocol
        const protocol = await this.syncService.getEntity<Protocol>(this.entityType, id);
        
        if (!protocol) {
          throw new Error(`Protocol with ID ${id} not found`);
        }
        
        // Convert to JSON and create a blob
        const json = JSON.stringify(protocol, null, 2);
        return new Blob([json], { type: 'application/json' });
      }
      
      // For other formats, we can't provide an offline alternative
      throw new Error(`Cannot export protocol as ${format} while offline`);
    }
  }
  
  /**
   * Imports a protocol from a file
   */
  async importProtocol(file: File) {
    try {
      // This operation requires server access for proper validation
      const formData = new FormData();
      formData.append('file', file);
      
      const response = await fetch(`${this.apiClient.baseURL}/${this.entityType}/import`, {
        method: 'POST',
        headers: {
          Authorization: `Bearer ${localStorage.getItem('auth_token') || ''}`,
        },
        body: formData,
      });
      
      if (!response.ok) {
        throw new Error(`Failed to import protocol: ${response.statusText}`);
      }
      
      const importedProtocol = await response.json();
      
      // Cache the imported protocol
      await this.syncService.saveEntity(this.entityType, importedProtocol.id, importedProtocol);
      
      return importedProtocol;
    } catch (error) {
      // For JSON files, we can provide an offline alternative
      if (file.type === 'application/json') {
        console.log('Importing protocol from JSON offline');
        
        try {
          // Read the file
          const text = await file.text();
          const data = JSON.parse(text);
          
          // Create a new protocol from the imported data
          const tempId = uuidv4() as UUID;
          const timestamp = new Date().toISOString();
          const deviceId = this.syncService.getDeviceId();
          
          // Create a valid protocol structure
          const importedProtocol: Protocol = {
            id: tempId,
            name: data.name || 'Imported Protocol',
            description: data.description || 'Imported while offline',
            version: data.version || '1.0.0',
            steps: Array.isArray(data.steps) ? data.steps : [],
            created_at: timestamp,
            updated_at: timestamp,
            created_by: deviceId,
            is_template: false,
            provenance: {
              created_by: deviceId,
              created_at: timestamp,
              source: 'offline-import',
              method: 'ui',
              version: '1.0',
              device_id: deviceId,
            },
            _temp_id: true, // Mark as temporary
          };
          
          // Ensure steps have proper IDs
          importedProtocol.steps = importedProtocol.steps.map((step, index) => ({
            ...step,
            id: step.id || uuidv4() as UUID,
            protocol_id: tempId,
            order: step.order || index + 1,
          }));
          
          // Save to offline storage
          await this.syncService.saveEntity(this.entityType, tempId, importedProtocol);
          
          return importedProtocol;
        } catch (parseError) {
          throw new Error(`Failed to parse JSON: ${parseError instanceof Error ? parseError.message : String(parseError)}`);
        }
      }
      
      // For other formats, we can't provide an offline alternative
      throw new Error('Cannot import protocol while offline except for JSON format');
    }
  }
  
  /**
   * Gets protocol templates
   */
  async getTemplates(params?: ProtocolListParams) {
    // Add is_template filter
    const templatesParams = {
      ...params,
      is_template: true,
    };
    
    // Use the standard getProtocols method with the is_template filter
    return this.getProtocols(templatesParams);
  }
  
  /**
   * Creates a new protocol from a template
   */
  async createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>) {
    try {
      // Try to create on the server first
      const response = await this.apiClient.post<Protocol>(
        `/${this.entityType}/create-from-template`,
        {
          template_id: templateId,
          customizations,
        }
      );
      
      // Cache the created protocol
      await this.syncService.saveEntity(this.entityType, response.id, response);
      
      return response;
    } catch (error) {
      // If API fails, create locally with a temporary ID
      console.log('Creating protocol from template offline');
      
      // Get the template
      const template = await this.syncService.getEntity<Protocol>(this.entityType, templateId);
      
      if (!template) {
        throw new Error(`Template with ID ${templateId} not found`);
      }
      
      // Create a new protocol based on the template
      const tempId = uuidv4() as UUID;
      const timestamp = new Date().toISOString();
      const deviceId = this.syncService.getDeviceId();
      
      // Create a new protocol from the template with customizations
      const newProtocol: Protocol = {
        ...template,
        ...customizations,
        id: tempId,
        is_template: false,
        version: '1.0.0',
        created_at: timestamp,
        updated_at: timestamp,
        created_by: customizations?.created_by || deviceId,
        provenance: {
          created_by: customizations?.created_by || deviceId,
          created_at: timestamp,
          source: 'offline-template',
          method: 'ui',
          version: '1.0',
          device_id: deviceId,
          references: [templateId.toString()],
        },
        _temp_id: true, // Mark as temporary
      };
      
      // Ensure steps have proper IDs
      newProtocol.steps = newProtocol.steps.map((step, index) => ({
        ...step,
        id: uuidv4() as UUID,
        protocol_id: tempId,
        order: step.order || index + 1,
      }));
      
      // Save to offline storage
      await this.syncService.saveEntity(this.entityType, tempId, newProtocol);
      
      return newProtocol;
    }
  }
  
  /**
   * Gets protocols from the library
   */
  async getLibraryProtocols(params?: ProtocolListParams) {
    try {
      // This operation requires server access
      const queryString = this.buildQueryParams(params);
      return await this.apiClient.get<{
        data: ProtocolLibraryItem[];
        total: number;
        page: number;
        per_page: number;
        total_pages: number;
      }>(`/${this.entityType}/library${queryString}`);
    } catch (error) {
      // For now, we can't provide a full offline library
      console.log('Using limited offline data for library protocols');
      
      // Get all local templates as a basic library
      const templatesParams = {
        ...params,
        is_template: true,
      };
      
      const templates = await this.getProtocols(templatesParams);
      
      // Convert to library items
      const libraryItems: ProtocolLibraryItem[] = templates.data.map(template => ({
        id: template.id,
        name: template.name,
        description: template.description,
        category: 'Offline',
        tags: template.tags || [],
        experiment_type: template.steps[0]?.parameters?.experiment_type || 'unknown',
        version: template.version,
        author: template.created_by,
        rating: 0,
        ratings_count: 0,
        used_count: 0,
        created_at: template.created_at,
        updated_at: template.updated_at,
      }));
      
      return {
        data: libraryItems,
        total: libraryItems.length,
        page: params?.page || 1,
        per_page: params?.per_page || 10,
        total_pages: Math.ceil(libraryItems.length / (params?.per_page || 10)),
      };
    }
  }
  
  /**
   * Gets a protocol from the library
   */
  async getLibraryProtocol(id: UUID) {
    try {
      // This operation requires server access
      const response = await this.apiClient.get<Protocol>(`/${this.entityType}/library/${id}`);
      
      // Cache the protocol
      await this.syncService.saveEntity(this.entityType, id, response);
      
      return response;
    } catch (error) {
      // If API fails, try to get from offline storage
      console.log(`Using offline data for library protocol ${id}`);
      
      const protocol = await this.syncService.getEntity<Protocol>(this.entityType, id);
      
      if (!protocol) {
        throw new Error(`Library protocol with ID ${id} not found`);
      }
      
      return protocol;
    }
  }
  
  /**
   * Publishes a protocol to the library
   */
  async publishToLibrary(protocolId: UUID, metadata: {
    category: string;
    tags: string[];
    description?: string;
  }) {
    try {
      // This operation requires server access
      const response = await this.apiClient.post<ProtocolLibraryItem>(
        `/${this.entityType}/${protocolId}/publish`,
        metadata
      );
      
      return response;
    } catch (error) {
      // This operation is not available offline
      throw new Error('Cannot publish protocol to library while offline');
    }
  }
  
  /**
   * Searches protocols
   */
  async searchProtocols(query: string, params?: ProtocolListParams) {
    try {
      // Try to search on the server first
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
      
      const response = await this.apiClient.get<{
        data: Protocol[];
        total: number;
      }>(`/${this.entityType}/search?${queryParams.toString()}`);
      
      // Cache each protocol
      for (const protocol of response.data) {
        await this.syncService.saveEntity(this.entityType, protocol.id, protocol);
      }
      
      return response;
    } catch (error) {
      // If API fails, search locally
      console.log('Searching protocols offline');
      
      // Get all protocols
      const protocols = await this.syncService.getAllEntities<Protocol>(this.entityType);
      
      // Perform basic search
      const searchQuery = query.toLowerCase();
      const results = protocols.filter(protocol => 
        protocol.name.toLowerCase().includes(searchQuery) ||
        (protocol.description && protocol.description.toLowerCase().includes(searchQuery)) ||
        protocol.steps.some(step => 
          step.name.toLowerCase().includes(searchQuery) ||
          (step.description && step.description.toLowerCase().includes(searchQuery))
        )
      );
      
      return {
        data: results,
        total: results.length,
      };
    }
  }
}