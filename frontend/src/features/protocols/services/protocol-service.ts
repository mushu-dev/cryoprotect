/**
 * Protocol service interface for managing experiment protocols
 */

import { UUID } from 'crypto';
import { Protocol, ProtocolStep, Provenance } from '../../experiments/services/experiment-service';

export interface ProtocolDiff {
  id: string;
  previous_version: string;
  current_version: string;
  changes: {
    added: ProtocolStep[];
    removed: ProtocolStep[];
    modified: {
      step_id: UUID;
      previous: Partial<ProtocolStep>;
      current: Partial<ProtocolStep>;
      changes: string[];
    }[];
  };
  parameter_changes?: {
    added: string[];
    removed: string[];
    modified: {
      name: string;
      previous: any;
      current: any;
    }[];
  };
  metadata_changes?: Record<string, {
    previous: any;
    current: any;
  }>;
  compatibility_analysis: {
    compatible: boolean;
    breaking_changes: string[];
    warnings: string[];
  };
  created_at: string;
}

export interface ValidationResult {
  valid: boolean;
  errors: {
    step_id?: UUID;
    field?: string;
    message: string;
    severity: 'error' | 'warning' | 'info';
  }[];
  warnings: {
    step_id?: UUID;
    field?: string;
    message: string;
  }[];
  suggestions: {
    step_id?: UUID;
    field?: string;
    message: string;
    suggested_value?: any;
  }[];
}

export interface ProtocolLibraryItem {
  id: UUID;
  name: string;
  description?: string;
  category: string;
  tags: string[];
  experiment_type: string;
  version: string;
  author: string;
  organization?: string;
  rating: number;
  ratings_count: number;
  used_count: number;
  thumbnail_url?: string;
  created_at: string;
  updated_at: string;
}

export interface ProtocolListParams {
  page?: number;
  per_page?: number;
  sort_by?: string;
  sort_order?: 'asc' | 'desc';
  experiment_type?: string;
  is_template?: boolean;
  author?: string;
  search?: string;
  tags?: string[];
  category?: string;
}

export interface ProtocolService {
  getProtocols(params?: ProtocolListParams): Promise<{
    data: Protocol[];
    total: number;
    page: number;
    per_page: number;
    total_pages: number;
  }>;
  
  getProtocol(id: UUID): Promise<Protocol>;
  
  getProtocolVersions(id: UUID): Promise<{
    protocol_id: UUID;
    versions: {
      version: string;
      id: UUID;
      created_at: string;
      created_by: string;
      description?: string;
    }[];
  }>;
  
  compareProtocolVersions(
    protocolId1: UUID,
    protocolId2: UUID
  ): Promise<ProtocolDiff>;
  
  createProtocol(protocol: Omit<Protocol, 'id' | 'provenance'>): Promise<Protocol>;
  
  updateProtocol(id: UUID, protocol: Partial<Protocol>): Promise<Protocol>;
  
  deleteProtocol(id: UUID): Promise<void>;
  
  createProtocolVersion(
    protocolId: UUID,
    changes: Partial<Protocol>,
    versionNotes?: string
  ): Promise<Protocol>;
  
  validateProtocol(protocol: Protocol): Promise<ValidationResult>;
  
  exportProtocol(
    id: UUID,
    format: 'json' | 'yaml' | 'pdf' | 'human-readable'
  ): Promise<Blob>;
  
  importProtocol(file: File): Promise<Protocol>;
  
  getTemplates(params?: ProtocolListParams): Promise<{
    data: Protocol[];
    total: number;
    page: number;
    per_page: number;
    total_pages: number;
  }>;
  
  createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>): Promise<Protocol>;
  
  getLibraryProtocols(params?: ProtocolListParams): Promise<{
    data: ProtocolLibraryItem[];
    total: number;
    page: number;
    per_page: number;
    total_pages: number;
  }>;
  
  getLibraryProtocol(id: UUID): Promise<Protocol>;
  
  publishToLibrary(protocolId: UUID, metadata: {
    category: string;
    tags: string[];
    description?: string;
  }): Promise<ProtocolLibraryItem>;
  
  searchProtocols(query: string, params?: ProtocolListParams): Promise<{
    data: Protocol[];
    total: number;
  }>;
}

export class ProtocolServiceImpl implements ProtocolService {
  private apiUrl: string;
  
  constructor(baseUrl: string = '/api') {
    this.apiUrl = `${baseUrl}/protocols`;
  }
  
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
    
    const response = await fetch(`${this.apiUrl}?${queryParams.toString()}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch protocols: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async getProtocol(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch protocol: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async getProtocolVersions(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}/versions`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch protocol versions: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async compareProtocolVersions(protocolId1: UUID, protocolId2: UUID) {
    const response = await fetch(`${this.apiUrl}/compare?id1=${protocolId1}&id2=${protocolId2}`);
    
    if (!response.ok) {
      throw new Error(`Failed to compare protocols: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async createProtocol(protocol: Omit<Protocol, 'id' | 'provenance'>) {
    const response = await fetch(this.apiUrl, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(protocol),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to create protocol: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async updateProtocol(id: UUID, protocol: Partial<Protocol>) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      method: 'PATCH',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(protocol),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to update protocol: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async deleteProtocol(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      method: 'DELETE',
    });
    
    if (!response.ok) {
      throw new Error(`Failed to delete protocol: ${response.statusText}`);
    }
  }
  
  async createProtocolVersion(protocolId: UUID, changes: Partial<Protocol>, versionNotes?: string) {
    const response = await fetch(`${this.apiUrl}/${protocolId}/versions`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        changes,
        version_notes: versionNotes,
      }),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to create protocol version: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async validateProtocol(protocol: Protocol) {
    const response = await fetch(`${this.apiUrl}/validate`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(protocol),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to validate protocol: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async exportProtocol(id: UUID, format: 'json' | 'yaml' | 'pdf' | 'human-readable') {
    const response = await fetch(`${this.apiUrl}/${id}/export?format=${format}`, {
      method: 'GET',
    });
    
    if (!response.ok) {
      throw new Error(`Failed to export protocol: ${response.statusText}`);
    }
    
    return await response.blob();
  }
  
  async importProtocol(file: File) {
    const formData = new FormData();
    formData.append('file', file);
    
    const response = await fetch(`${this.apiUrl}/import`, {
      method: 'POST',
      body: formData,
    });
    
    if (!response.ok) {
      throw new Error(`Failed to import protocol: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
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
    
    const response = await fetch(`${this.apiUrl}?${queryParams.toString()}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch protocol templates: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>) {
    const response = await fetch(`${this.apiUrl}/create-from-template`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        template_id: templateId,
        customizations,
      }),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to create protocol from template: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
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
    
    const response = await fetch(`${this.apiUrl}/library?${queryParams.toString()}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch library protocols: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async getLibraryProtocol(id: UUID) {
    const response = await fetch(`${this.apiUrl}/library/${id}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch library protocol: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async publishToLibrary(protocolId: UUID, metadata: {
    category: string;
    tags: string[];
    description?: string;
  }) {
    const response = await fetch(`${this.apiUrl}/${protocolId}/publish`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(metadata),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to publish protocol to library: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
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
    
    const response = await fetch(`${this.apiUrl}/search?${queryParams.toString()}`);
    
    if (!response.ok) {
      throw new Error(`Failed to search protocols: ${response.statusText}`);
    }
    
    return await response.json();
  }
}

export class MockProtocolService implements ProtocolService {
  private protocols: Protocol[] = [];
  private libraryProtocols: ProtocolLibraryItem[] = [];
  
  constructor() {
    // Initialize with mock data if needed
  }
  
  async getProtocols(params?: ProtocolListParams) {
    const page = params?.page || 1;
    const perPage = params?.per_page || 10;
    const startIndex = (page - 1) * perPage;
    const endIndex = startIndex + perPage;
    
    let filteredProtocols = [...this.protocols];
    
    // Apply filters
    if (params) {
      if (params.experiment_type) {
        filteredProtocols = filteredProtocols.filter(p => p.steps.some(s => 
          s.parameters && s.parameters.experiment_type === params.experiment_type
        ));
      }
      
      if (params.is_template !== undefined) {
        filteredProtocols = filteredProtocols.filter(p => p.is_template === params.is_template);
      }
      
      if (params.author) {
        filteredProtocols = filteredProtocols.filter(p => p.created_by === params.author);
      }
      
      if (params.tags && params.tags.length > 0) {
        filteredProtocols = filteredProtocols.filter(p => 
          p.tags && params.tags?.some(tag => p.tags?.includes(tag))
        );
      }
      
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
    }
    
    // Apply sorting
    if (params?.sort_by) {
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
    
    const paginatedProtocols = filteredProtocols.slice(startIndex, endIndex);
    
    return {
      data: paginatedProtocols,
      total: filteredProtocols.length,
      page,
      per_page: perPage,
      total_pages: Math.ceil(filteredProtocols.length / perPage),
    };
  }
  
  async getProtocol(id: UUID) {
    const protocol = this.protocols.find(p => p.id === id);
    
    if (!protocol) {
      throw new Error(`Protocol with ID ${id} not found`);
    }
    
    return { ...protocol };
  }
  
  async getProtocolVersions(id: UUID) {
    const versions = this.protocols
      .filter(p => p.id === id || p.parent_version_id === id)
      .map(p => ({
        version: p.version,
        id: p.id,
        created_at: p.created_at,
        created_by: p.created_by,
        description: p.description,
      }));
    
    return {
      protocol_id: id,
      versions,
    };
  }
  
  async compareProtocolVersions() {
    // Return mock comparison result
    return {
      id: crypto.randomUUID(),
      previous_version: '1.0.0',
      current_version: '1.1.0',
      changes: {
        added: [],
        removed: [],
        modified: [
          {
            step_id: crypto.randomUUID() as UUID,
            previous: {
              temperature: 20,
              temperature_unit: 'C',
            },
            current: {
              temperature: 25,
              temperature_unit: 'C',
            },
            changes: ['temperature changed from 20°C to 25°C'],
          },
        ],
      },
      parameter_changes: {
        added: ['new_parameter'],
        removed: [],
        modified: [
          {
            name: 'existing_parameter',
            previous: 'old value',
            current: 'new value',
          },
        ],
      },
      compatibility_analysis: {
        compatible: true,
        breaking_changes: [],
        warnings: ['Temperature increase may affect outcome'],
      },
      created_at: new Date().toISOString(),
    };
  }
  
  async createProtocol(protocol: Omit<Protocol, 'id' | 'provenance'>) {
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    const newProtocol: Protocol = {
      ...protocol as any,
      id,
      version: protocol.version || '1.0.0',
      created_at: timestamp,
      updated_at: timestamp,
      provenance: {
        created_by: protocol.created_by,
        created_at: timestamp,
        source: 'manual-entry',
        method: 'ui',
        version: '1.0',
      },
    };
    
    this.protocols.push(newProtocol);
    
    return newProtocol;
  }
  
  async updateProtocol(id: UUID, protocol: Partial<Protocol>) {
    const index = this.protocols.findIndex(p => p.id === id);
    
    if (index === -1) {
      throw new Error(`Protocol with ID ${id} not found`);
    }
    
    const updatedProtocol = {
      ...this.protocols[index],
      ...protocol,
      updated_at: new Date().toISOString(),
    };
    
    this.protocols[index] = updatedProtocol;
    
    return updatedProtocol;
  }
  
  async deleteProtocol(id: UUID) {
    const index = this.protocols.findIndex(p => p.id === id);
    
    if (index === -1) {
      throw new Error(`Protocol with ID ${id} not found`);
    }
    
    this.protocols.splice(index, 1);
  }
  
  async createProtocolVersion(protocolId: UUID, changes: Partial<Protocol>) {
    const protocol = this.protocols.find(p => p.id === protocolId);
    
    if (!protocol) {
      throw new Error(`Protocol with ID ${protocolId} not found`);
    }
    
    // Parse current version and increment minor version
    const versionParts = protocol.version.split('.');
    const newVersion = `${versionParts[0]}.${parseInt(versionParts[1]) + 1}.0`;
    
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    const newVersionProtocol: Protocol = {
      ...protocol,
      ...changes,
      id,
      version: newVersion,
      parent_version_id: protocolId,
      created_at: timestamp,
      updated_at: timestamp,
      provenance: {
        created_by: protocol.created_by,
        created_at: timestamp,
        source: 'version-creation',
        method: 'ui',
        version: '1.0',
      },
    };
    
    this.protocols.push(newVersionProtocol);
    
    return newVersionProtocol;
  }
  
  async validateProtocol() {
    // Return mock validation result
    return {
      valid: true,
      errors: [],
      warnings: [
        {
          step_id: crypto.randomUUID() as UUID,
          field: 'temperature',
          message: 'Temperature is outside the recommended range for this experiment type',
        },
      ],
      suggestions: [
        {
          step_id: crypto.randomUUID() as UUID,
          field: 'duration',
          message: 'Consider increasing duration for better results',
          suggested_value: 30,
        },
      ],
    };
  }
  
  async exportProtocol() {
    // Return mock blob
    return new Blob(['Mock protocol export data'], { type: 'application/json' });
  }
  
  async importProtocol() {
    // Return mock imported protocol
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    const importedProtocol: Protocol = {
      id,
      name: 'Imported Protocol',
      description: 'This is an imported protocol',
      version: '1.0.0',
      steps: [
        {
          id: crypto.randomUUID() as UUID,
          protocol_id: id,
          name: 'Step 1',
          description: 'First step of imported protocol',
          order: 1,
          duration: 10,
          duration_unit: 'minutes',
          temperature: 20,
          temperature_unit: 'C',
        },
      ],
      created_at: timestamp,
      updated_at: timestamp,
      created_by: 'import-user',
      is_template: false,
      provenance: {
        created_by: 'import-user',
        created_at: timestamp,
        source: 'file-import',
        method: 'ui',
        version: '1.0',
      },
    };
    
    this.protocols.push(importedProtocol);
    
    return importedProtocol;
  }
  
  async getTemplates(params?: ProtocolListParams) {
    const templatesParams = {
      ...params,
      is_template: true,
    };
    
    return this.getProtocols(templatesParams);
  }
  
  async createFromTemplate(templateId: UUID, customizations?: Partial<Protocol>) {
    const template = this.protocols.find(p => p.id === templateId);
    
    if (!template) {
      throw new Error(`Template with ID ${templateId} not found`);
    }
    
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    // Assign new IDs to steps to avoid conflicts
    const steps = template.steps.map(step => ({
      ...step,
      id: crypto.randomUUID() as UUID,
      protocol_id: id,
    }));
    
    const newProtocol: Protocol = {
      ...template,
      ...customizations,
      id,
      steps,
      is_template: false,
      created_at: timestamp,
      updated_at: timestamp,
      version: '1.0.0',
      provenance: {
        created_by: customizations?.created_by || 'unknown',
        created_at: timestamp,
        source: 'template',
        method: 'ui',
        version: '1.0',
        references: [templateId.toString()],
      },
    };
    
    this.protocols.push(newProtocol);
    
    return newProtocol;
  }
  
  async getLibraryProtocols(params?: ProtocolListParams) {
    const page = params?.page || 1;
    const perPage = params?.per_page || 10;
    const startIndex = (page - 1) * perPage;
    const endIndex = startIndex + perPage;
    
    let filteredProtocols = [...this.libraryProtocols];
    
    // Apply filters
    if (params) {
      if (params.experiment_type) {
        filteredProtocols = filteredProtocols.filter(p => p.experiment_type === params.experiment_type);
      }
      
      if (params.category) {
        filteredProtocols = filteredProtocols.filter(p => p.category === params.category);
      }
      
      if (params.author) {
        filteredProtocols = filteredProtocols.filter(p => p.author === params.author);
      }
      
      if (params.tags && params.tags.length > 0) {
        filteredProtocols = filteredProtocols.filter(p => 
          params.tags?.some(tag => p.tags.includes(tag))
        );
      }
      
      if (params.search) {
        const search = params.search.toLowerCase();
        filteredProtocols = filteredProtocols.filter(p => 
          p.name.toLowerCase().includes(search) || 
          (p.description && p.description.toLowerCase().includes(search))
        );
      }
    }
    
    // Apply sorting
    if (params?.sort_by) {
      const sortBy = params.sort_by as keyof ProtocolLibraryItem;
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
    
    const paginatedProtocols = filteredProtocols.slice(startIndex, endIndex);
    
    return {
      data: paginatedProtocols,
      total: filteredProtocols.length,
      page,
      per_page: perPage,
      total_pages: Math.ceil(filteredProtocols.length / perPage),
    };
  }
  
  async getLibraryProtocol(id: UUID) {
    const libraryItem = this.libraryProtocols.find(p => p.id === id);
    
    if (!libraryItem) {
      throw new Error(`Library protocol with ID ${id} not found`);
    }
    
    // Mock the full protocol data
    const protocol: Protocol = {
      id,
      name: libraryItem.name,
      description: libraryItem.description,
      version: libraryItem.version,
      steps: [],
      created_at: libraryItem.created_at,
      updated_at: libraryItem.updated_at,
      created_by: libraryItem.author,
      is_template: true,
      tags: libraryItem.tags,
      provenance: {
        created_by: libraryItem.author,
        created_at: libraryItem.created_at,
        source: 'library',
        method: 'ui',
        version: '1.0',
      },
    };
    
    return protocol;
  }
  
  async publishToLibrary(protocolId: UUID, metadata: {
    category: string;
    tags: string[];
    description?: string;
  }) {
    const protocol = this.protocols.find(p => p.id === protocolId);
    
    if (!protocol) {
      throw new Error(`Protocol with ID ${protocolId} not found`);
    }
    
    const libraryItem: ProtocolLibraryItem = {
      id: protocolId,
      name: protocol.name,
      description: metadata.description || protocol.description,
      category: metadata.category,
      tags: metadata.tags,
      experiment_type: protocol.steps.length > 0 && protocol.steps[0].parameters?.experiment_type || 'unknown',
      version: protocol.version,
      author: protocol.created_by,
      rating: 0,
      ratings_count: 0,
      used_count: 0,
      created_at: protocol.created_at,
      updated_at: new Date().toISOString(),
    };
    
    this.libraryProtocols.push(libraryItem);
    
    return libraryItem;
  }
  
  async searchProtocols(query: string, params?: ProtocolListParams) {
    const searchParams = {
      ...params,
      search: query,
    };
    
    const result = await this.getProtocols(searchParams);
    
    return {
      data: result.data,
      total: result.total,
    };
  }
}