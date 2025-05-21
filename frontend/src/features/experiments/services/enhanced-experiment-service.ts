/**
 * Enhanced Experiment Service with resilient API layer integration
 */

import { UUID } from 'crypto';
import { 
  ExperimentService, 
  Experiment, 
  ExperimentResult, 
  TimeSeries,
  ExperimentListParams,
  ExperimentAnalysisResult
} from './experiment-service';
import { experimentalDataApi, getApiService } from '../../api/api-service';

/**
 * Enhanced experiment service implementation that uses the resilient API layer
 * with circuit breaker and fallback mechanisms
 */
export class ResilienceEnhancedExperimentService implements ExperimentService {
  constructor(private fallbackService?: ExperimentService) {}

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
  private convertParams(params?: ExperimentListParams): Record<string, any> {
    if (!params) return {};
    
    const apiParams: Record<string, any> = {};
    
    for (const [key, value] of Object.entries(params)) {
      if (value !== undefined) {
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
   * Get experiments with pagination and filtering
   */
  async getExperiments(params?: ExperimentListParams) {
    try {
      const apiParams = this.convertParams(params);
      const response = await this.getApi().getExperiments(apiParams);
      
      return {
        data: response.data || [],
        total: response.total || 0,
        page: response.page || 1,
        per_page: response.per_page || 10,
        total_pages: response.total_pages || 1
      };
    } catch (error) {
      console.error('Error fetching experiments:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getExperiments');
        return this.fallbackService.getExperiments(params);
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
   * Get experiment by ID
   */
  async getExperiment(id: UUID) {
    try {
      return await this.getApi().getExperiment(id);
    } catch (error) {
      console.error(`Error fetching experiment ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getExperiment');
        return this.fallbackService.getExperiment(id);
      }
      
      throw error;
    }
  }

  /**
   * Create new experiment
   */
  async createExperiment(experiment: Omit<Experiment, 'id' | 'provenance'>) {
    try {
      return await this.getApi().createExperiment(experiment);
    } catch (error) {
      console.error('Error creating experiment:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for createExperiment');
        return this.fallbackService.createExperiment(experiment);
      }
      
      throw error;
    }
  }

  /**
   * Update existing experiment
   */
  async updateExperiment(id: UUID, experiment: Partial<Experiment>) {
    try {
      return await this.getApi().updateExperiment(id, experiment);
    } catch (error) {
      console.error(`Error updating experiment ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for updateExperiment');
        return this.fallbackService.updateExperiment(id, experiment);
      }
      
      throw error;
    }
  }

  /**
   * Delete experiment
   */
  async deleteExperiment(id: UUID) {
    try {
      await this.getApi().deleteExperiment(id);
    } catch (error) {
      console.error(`Error deleting experiment ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for deleteExperiment');
        return this.fallbackService.deleteExperiment(id);
      }
      
      throw error;
    }
  }

  /**
   * Get experiment results
   */
  async getExperimentResults(experimentId: UUID) {
    try {
      return await this.getApi().getExperimentResults(experimentId);
    } catch (error) {
      console.error(`Error fetching results for experiment ${experimentId}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getExperimentResults');
        return this.fallbackService.getExperimentResults(experimentId);
      }
      
      return [];
    }
  }

  /**
   * Add experiment result
   */
  async addExperimentResult(result: Omit<ExperimentResult, 'id' | 'provenance'>) {
    try {
      return await this.getApi().addExperimentResult(result);
    } catch (error) {
      console.error('Error adding experiment result:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for addExperimentResult');
        return this.fallbackService.addExperimentResult(result);
      }
      
      throw error;
    }
  }

  /**
   * Update experiment result
   */
  async updateExperimentResult(id: UUID, result: Partial<ExperimentResult>) {
    try {
      return await this.getApi().updateExperimentResult(id, result);
    } catch (error) {
      console.error(`Error updating experiment result ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for updateExperimentResult');
        return this.fallbackService.updateExperimentResult(id, result);
      }
      
      throw error;
    }
  }

  /**
   * Delete experiment result
   */
  async deleteExperimentResult(id: UUID) {
    try {
      await this.getApi().deleteExperimentResult(id);
    } catch (error) {
      console.error(`Error deleting experiment result ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for deleteExperimentResult');
        return this.fallbackService.deleteExperimentResult(id);
      }
      
      throw error;
    }
  }

  /**
   * Get experiment time series data
   */
  async getExperimentTimeSeries(experimentId: UUID, parameter?: string) {
    try {
      return await this.getApi().getExperimentTimeSeries(experimentId, parameter);
    } catch (error) {
      console.error(`Error fetching time series for experiment ${experimentId}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for getExperimentTimeSeries');
        return this.fallbackService.getExperimentTimeSeries(experimentId, parameter);
      }
      
      return [];
    }
  }

  /**
   * Analyze experiments with statistical methods
   */
  async analyzeExperiments(
    experimentIds: UUID[],
    parameters?: {
      compare_with?: UUID[];
      analysis_type?: string[];
      time_range?: [string, string];
    }
  ) {
    try {
      return await this.getApi().analyzeExperiments(experimentIds, parameters);
    } catch (error) {
      console.error('Error analyzing experiments:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for analyzeExperiments');
        return this.fallbackService.analyzeExperiments(experimentIds, parameters);
      }
      
      // Return basic empty analysis if no fallback
      return {
        summary: {
          total: 0,
          successful: 0,
          failed: 0,
          success_rate: 0,
        },
        statistics: {
          viability: {
            mean: 0,
            median: 0,
            std_dev: 0,
            min: 0,
            max: 0,
            quartiles: [0, 0, 0],
          },
          recovery: {
            mean: 0,
            median: 0,
            std_dev: 0,
            min: 0,
            max: 0,
            quartiles: [0, 0, 0],
          },
        },
        trends: [],
      } as ExperimentAnalysisResult;
    }
  }

  /**
   * Export experiment data in various formats
   */
  async exportExperiment(
    id: UUID,
    format: 'json' | 'csv' | 'excel' | 'pdf'
  ) {
    try {
      return await this.getApi().exportExperiment(id, format);
    } catch (error) {
      console.error(`Error exporting experiment ${id}:`, error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for exportExperiment');
        return this.fallbackService.exportExperiment(id, format);
      }
      
      // Return empty blob if no fallback
      return new Blob([`Error exporting experiment ${id}`], { type: 'text/plain' });
    }
  }

  /**
   * Import experiment from file
   */
  async importExperiment(file: File) {
    try {
      return await this.getApi().importExperiment(file);
    } catch (error) {
      console.error('Error importing experiment:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for importExperiment');
        return this.fallbackService.importExperiment(file);
      }
      
      throw error;
    }
  }

  /**
   * Search experiments by query string
   */
  async searchExperiments(query: string, params?: ExperimentListParams) {
    try {
      return await this.getApi().searchExperiments(query, params);
    } catch (error) {
      console.error('Error searching experiments:', error);
      
      // Fallback to mock service if available
      if (this.fallbackService) {
        console.warn('Using fallback service for searchExperiments');
        return this.fallbackService.searchExperiments(query, params);
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
 * Create a resilient experiment service with automatic fallback
 */
export const createResilienceEnhancedExperimentService = (fallbackService?: ExperimentService) => {
  return new ResilienceEnhancedExperimentService(fallbackService);
};