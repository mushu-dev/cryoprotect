/**
 * Experiment service interface for managing experiments and their results
 */

import { UUID } from 'crypto';

export interface Uncertainty {
  value: number;
  type: 'standard' | 'expanded' | 'range';
  confidence: number;
  distribution: 'normal' | 'uniform' | 'triangular';
}

export interface Provenance {
  created_by: string;
  created_at: string;
  source: string;
  method: string;
  version: string;
  references?: string[];
}

export interface ExperimentResult {
  id: UUID;
  experiment_id: UUID;
  tissue_type_id: UUID;
  molecule_id?: UUID;
  mixture_id?: UUID;
  concentration?: number;
  concentration_unit?: string;
  viability_percentage?: number;
  recovery_rate?: number;
  functionality_score?: number;
  uncertainty?: {
    viability_percentage?: Uncertainty;
    recovery_rate?: Uncertainty;
    functionality_score?: Uncertainty;
    [key: string]: Uncertainty | undefined;
  };
  result_details?: Record<string, any>;
  notes?: string;
  protocol_step_id?: UUID;
  timestamp: string;
  provenance: Provenance;
}

export interface TimeSeries {
  id: UUID;
  experiment_id: UUID;
  result_id?: UUID;
  parameter: string;
  unit: string;
  data_points: {
    time: number;
    value: number;
    uncertainty?: Uncertainty;
  }[];
  start_time: string;
  end_time: string;
  notes?: string;
  provenance: Provenance;
}

export interface Protocol {
  id: UUID;
  name: string;
  description?: string;
  version: string;
  parent_version_id?: UUID;
  steps: ProtocolStep[];
  parameters?: Record<string, any>;
  tags?: string[];
  created_at: string;
  updated_at: string;
  created_by: string;
  is_template: boolean;
  provenance: Provenance;
}

export interface ProtocolStep {
  id: UUID;
  protocol_id: UUID;
  name: string;
  description?: string;
  order: number;
  duration?: number;
  duration_unit?: string;
  temperature?: number;
  temperature_unit?: string;
  parameters?: Record<string, any>;
  equipment?: string[];
  alerts?: {
    condition: string;
    message: string;
    severity: 'info' | 'warning' | 'critical';
  }[];
}

export interface Experiment {
  id: UUID;
  name: string;
  description?: string;
  protocol_id: UUID;
  protocol_version?: string;
  tissue_type_id: UUID;
  experiment_type: string;
  start_date: string;
  end_date?: string;
  status: 'planned' | 'in_progress' | 'completed' | 'aborted' | 'failed';
  researcher: string;
  lab_id?: string;
  equipment?: string[];
  environmental_conditions?: Record<string, any>;
  notes?: string;
  tags?: string[];
  results?: ExperimentResult[];
  time_series?: TimeSeries[];
  protocol?: Protocol;
  provenance: Provenance;
}

export interface ExperimentListParams {
  page?: number;
  per_page?: number;
  sort_by?: string;
  sort_order?: 'asc' | 'desc';
  status?: string;
  experiment_type?: string;
  researcher?: string;
  tissue_type_id?: UUID;
  date_from?: string;
  date_to?: string;
  tags?: string[];
  search?: string;
}

export interface ExperimentAnalysisResult {
  summary: {
    total: number;
    successful: number;
    failed: number;
    success_rate: number;
  };
  statistics: {
    viability: {
      mean: number;
      median: number;
      std_dev: number;
      min: number;
      max: number;
      quartiles: [number, number, number];
    };
    recovery: {
      mean: number;
      median: number;
      std_dev: number;
      min: number;
      max: number;
      quartiles: [number, number, number];
    };
  };
  trends: {
    parameter: string;
    timestamps: string[];
    values: number[];
    trend_line?: number[];
  }[];
  comparisons?: Record<string, any>;
}

export interface ExperimentService {
  getExperiments(params?: ExperimentListParams): Promise<{
    data: Experiment[];
    total: number;
    page: number;
    per_page: number;
    total_pages: number;
  }>;
  
  getExperiment(id: UUID): Promise<Experiment>;
  
  createExperiment(experiment: Omit<Experiment, 'id' | 'provenance'>): Promise<Experiment>;
  
  updateExperiment(id: UUID, experiment: Partial<Experiment>): Promise<Experiment>;
  
  deleteExperiment(id: UUID): Promise<void>;
  
  getExperimentResults(experimentId: UUID): Promise<ExperimentResult[]>;
  
  addExperimentResult(result: Omit<ExperimentResult, 'id' | 'provenance'>): Promise<ExperimentResult>;
  
  updateExperimentResult(id: UUID, result: Partial<ExperimentResult>): Promise<ExperimentResult>;
  
  deleteExperimentResult(id: UUID): Promise<void>;
  
  getExperimentTimeSeries(experimentId: UUID, parameter?: string): Promise<TimeSeries[]>;
  
  analyzeExperiments(
    experimentIds: UUID[],
    parameters?: {
      compare_with?: UUID[];
      analysis_type?: string[];
      time_range?: [string, string];
    }
  ): Promise<ExperimentAnalysisResult>;
  
  exportExperiment(
    id: UUID,
    format: 'json' | 'csv' | 'excel' | 'pdf'
  ): Promise<Blob>;
  
  importExperiment(file: File): Promise<Experiment>;

  searchExperiments(query: string, params?: ExperimentListParams): Promise<{
    data: Experiment[];
    total: number;
  }>;
}

export class ExperimentServiceImpl implements ExperimentService {
  private apiUrl: string;
  
  constructor(baseUrl: string = '/api') {
    this.apiUrl = `${baseUrl}/experiments`;
  }
  
  async getExperiments(params?: ExperimentListParams) {
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
      throw new Error(`Failed to fetch experiments: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async getExperiment(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch experiment: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async createExperiment(experiment: Omit<Experiment, 'id' | 'provenance'>) {
    const response = await fetch(this.apiUrl, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(experiment),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to create experiment: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async updateExperiment(id: UUID, experiment: Partial<Experiment>) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      method: 'PATCH',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(experiment),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to update experiment: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async deleteExperiment(id: UUID) {
    const response = await fetch(`${this.apiUrl}/${id}`, {
      method: 'DELETE',
    });
    
    if (!response.ok) {
      throw new Error(`Failed to delete experiment: ${response.statusText}`);
    }
  }
  
  async getExperimentResults(experimentId: UUID) {
    const response = await fetch(`${this.apiUrl}/${experimentId}/results`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch experiment results: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async addExperimentResult(result: Omit<ExperimentResult, 'id' | 'provenance'>) {
    const response = await fetch(`${this.apiUrl}/${result.experiment_id}/results`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(result),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to add experiment result: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async updateExperimentResult(id: UUID, result: Partial<ExperimentResult>) {
    const response = await fetch(`${this.apiUrl}/results/${id}`, {
      method: 'PATCH',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(result),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to update experiment result: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async deleteExperimentResult(id: UUID) {
    const response = await fetch(`${this.apiUrl}/results/${id}`, {
      method: 'DELETE',
    });
    
    if (!response.ok) {
      throw new Error(`Failed to delete experiment result: ${response.statusText}`);
    }
  }
  
  async getExperimentTimeSeries(experimentId: UUID, parameter?: string) {
    const queryParams = new URLSearchParams();
    
    if (parameter) {
      queryParams.append('parameter', parameter);
    }
    
    const response = await fetch(`${this.apiUrl}/${experimentId}/timeseries?${queryParams.toString()}`);
    
    if (!response.ok) {
      throw new Error(`Failed to fetch experiment time series: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async analyzeExperiments(
    experimentIds: UUID[],
    parameters?: {
      compare_with?: UUID[];
      analysis_type?: string[];
      time_range?: [string, string];
    }
  ) {
    const response = await fetch(`${this.apiUrl}/analyze`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        experiment_ids: experimentIds,
        ...parameters,
      }),
    });
    
    if (!response.ok) {
      throw new Error(`Failed to analyze experiments: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async exportExperiment(id: UUID, format: 'json' | 'csv' | 'excel' | 'pdf') {
    const response = await fetch(`${this.apiUrl}/${id}/export?format=${format}`, {
      method: 'GET',
    });
    
    if (!response.ok) {
      throw new Error(`Failed to export experiment: ${response.statusText}`);
    }
    
    return await response.blob();
  }
  
  async importExperiment(file: File) {
    const formData = new FormData();
    formData.append('file', file);
    
    const response = await fetch(`${this.apiUrl}/import`, {
      method: 'POST',
      body: formData,
    });
    
    if (!response.ok) {
      throw new Error(`Failed to import experiment: ${response.statusText}`);
    }
    
    return await response.json();
  }
  
  async searchExperiments(query: string, params?: ExperimentListParams) {
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
      throw new Error(`Failed to search experiments: ${response.statusText}`);
    }
    
    return await response.json();
  }
}

export class MockExperimentService implements ExperimentService {
  private experiments: Experiment[] = [];
  private results: ExperimentResult[] = [];
  private timeSeries: TimeSeries[] = [];
  
  constructor() {
    // Initialize with mock data if needed
  }
  
  async getExperiments(params?: ExperimentListParams) {
    const page = params?.page || 1;
    const perPage = params?.per_page || 10;
    const startIndex = (page - 1) * perPage;
    const endIndex = startIndex + perPage;
    
    let filteredExperiments = [...this.experiments];
    
    // Apply filters
    if (params) {
      if (params.status) {
        filteredExperiments = filteredExperiments.filter(e => e.status === params.status);
      }
      
      if (params.experiment_type) {
        filteredExperiments = filteredExperiments.filter(e => e.experiment_type === params.experiment_type);
      }
      
      if (params.researcher) {
        filteredExperiments = filteredExperiments.filter(e => e.researcher === params.researcher);
      }
      
      if (params.tissue_type_id) {
        filteredExperiments = filteredExperiments.filter(e => e.tissue_type_id === params.tissue_type_id);
      }
      
      if (params.date_from) {
        const fromDate = new Date(params.date_from);
        filteredExperiments = filteredExperiments.filter(e => new Date(e.start_date) >= fromDate);
      }
      
      if (params.date_to) {
        const toDate = new Date(params.date_to);
        filteredExperiments = filteredExperiments.filter(e => new Date(e.start_date) <= toDate);
      }
      
      if (params.tags && params.tags.length > 0) {
        filteredExperiments = filteredExperiments.filter(e => 
          e.tags && params.tags?.some(tag => e.tags?.includes(tag))
        );
      }
      
      if (params.search) {
        const search = params.search.toLowerCase();
        filteredExperiments = filteredExperiments.filter(e => 
          e.name.toLowerCase().includes(search) || 
          (e.description && e.description.toLowerCase().includes(search))
        );
      }
    }
    
    // Apply sorting
    if (params?.sort_by) {
      const sortBy = params.sort_by as keyof Experiment;
      const sortOrder = params.sort_order === 'desc' ? -1 : 1;
      
      filteredExperiments.sort((a, b) => {
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
    
    const paginatedExperiments = filteredExperiments.slice(startIndex, endIndex);
    
    return {
      data: paginatedExperiments,
      total: filteredExperiments.length,
      page,
      per_page: perPage,
      total_pages: Math.ceil(filteredExperiments.length / perPage),
    };
  }
  
  async getExperiment(id: UUID) {
    const experiment = this.experiments.find(e => e.id === id);
    
    if (!experiment) {
      throw new Error(`Experiment with ID ${id} not found`);
    }
    
    return { ...experiment };
  }
  
  async createExperiment(experiment: Omit<Experiment, 'id' | 'provenance'>) {
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    const newExperiment: Experiment = {
      ...experiment as any,
      id,
      provenance: {
        created_by: 'mock-user',
        created_at: timestamp,
        source: 'manual-entry',
        method: 'ui',
        version: '1.0',
      },
    };
    
    this.experiments.push(newExperiment);
    
    return newExperiment;
  }
  
  async updateExperiment(id: UUID, experiment: Partial<Experiment>) {
    const index = this.experiments.findIndex(e => e.id === id);
    
    if (index === -1) {
      throw new Error(`Experiment with ID ${id} not found`);
    }
    
    const updatedExperiment = {
      ...this.experiments[index],
      ...experiment,
      updated_at: new Date().toISOString(),
    };
    
    this.experiments[index] = updatedExperiment;
    
    return updatedExperiment;
  }
  
  async deleteExperiment(id: UUID) {
    const index = this.experiments.findIndex(e => e.id === id);
    
    if (index === -1) {
      throw new Error(`Experiment with ID ${id} not found`);
    }
    
    this.experiments.splice(index, 1);
    
    // Also remove related results and time series
    this.results = this.results.filter(r => r.experiment_id !== id);
    this.timeSeries = this.timeSeries.filter(ts => ts.experiment_id !== id);
  }
  
  async getExperimentResults(experimentId: UUID) {
    return this.results.filter(r => r.experiment_id === experimentId);
  }
  
  async addExperimentResult(result: Omit<ExperimentResult, 'id' | 'provenance'>) {
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    const newResult: ExperimentResult = {
      ...result as any,
      id,
      provenance: {
        created_by: 'mock-user',
        created_at: timestamp,
        source: 'manual-entry',
        method: 'ui',
        version: '1.0',
      },
      timestamp: timestamp,
    };
    
    this.results.push(newResult);
    
    return newResult;
  }
  
  async updateExperimentResult(id: UUID, result: Partial<ExperimentResult>) {
    const index = this.results.findIndex(r => r.id === id);
    
    if (index === -1) {
      throw new Error(`Experiment result with ID ${id} not found`);
    }
    
    const updatedResult = {
      ...this.results[index],
      ...result,
    };
    
    this.results[index] = updatedResult;
    
    return updatedResult;
  }
  
  async deleteExperimentResult(id: UUID) {
    const index = this.results.findIndex(r => r.id === id);
    
    if (index === -1) {
      throw new Error(`Experiment result with ID ${id} not found`);
    }
    
    this.results.splice(index, 1);
  }
  
  async getExperimentTimeSeries(experimentId: UUID, parameter?: string) {
    let filteredTimeSeries = this.timeSeries.filter(ts => ts.experiment_id === experimentId);
    
    if (parameter) {
      filteredTimeSeries = filteredTimeSeries.filter(ts => ts.parameter === parameter);
    }
    
    return filteredTimeSeries;
  }
  
  async analyzeExperiments() {
    // Return mock analysis result
    return {
      summary: {
        total: 10,
        successful: 8,
        failed: 2,
        success_rate: 0.8,
      },
      statistics: {
        viability: {
          mean: 75.2,
          median: 78.5,
          std_dev: 12.3,
          min: 45.6,
          max: 91.2,
          quartiles: [65.4, 78.5, 85.6],
        },
        recovery: {
          mean: 68.7,
          median: 71.2,
          std_dev: 15.6,
          min: 38.9,
          max: 88.4,
          quartiles: [58.2, 71.2, 82.3],
        },
      },
      trends: [
        {
          parameter: 'viability',
          timestamps: [
            '2023-01-01T00:00:00Z',
            '2023-02-01T00:00:00Z',
            '2023-03-01T00:00:00Z',
            '2023-04-01T00:00:00Z',
          ],
          values: [65.2, 72.3, 76.8, 80.1],
          trend_line: [64.5, 70.2, 75.9, 81.6],
        },
      ],
    };
  }
  
  async exportExperiment() {
    // Return mock blob
    return new Blob(['Mock export data'], { type: 'application/json' });
  }
  
  async importExperiment() {
    // Return mock imported experiment
    const id = crypto.randomUUID() as UUID;
    const timestamp = new Date().toISOString();
    
    const importedExperiment: Experiment = {
      id,
      name: 'Imported Experiment',
      description: 'This is an imported experiment',
      protocol_id: crypto.randomUUID() as UUID,
      tissue_type_id: crypto.randomUUID() as UUID,
      experiment_type: 'vitrification',
      start_date: timestamp,
      status: 'completed',
      researcher: 'Imported Researcher',
      provenance: {
        created_by: 'import-user',
        created_at: timestamp,
        source: 'file-import',
        method: 'ui',
        version: '1.0',
      },
    };
    
    this.experiments.push(importedExperiment);
    
    return importedExperiment;
  }
  
  async searchExperiments(query: string, params?: ExperimentListParams) {
    const searchParams = {
      ...params,
      search: query,
    };
    
    const result = await this.getExperiments(searchParams);
    
    return {
      data: result.data,
      total: result.total,
    };
  }
}