/**
 * Protocol Integration Service for managing protocol-experiment relationships
 */

import { UUID } from 'crypto';
import { Protocol } from '../../experiments/services/experiment-service';
import { Id } from '../../../convex/_generated/dataModel';
import { ConvexClient } from '@convex/react';
import { api } from '../../../convex/_generated/api';

export interface ExperimentSummary {
  id: string;
  title: string;
  description?: string;
  status: string;
  datePerformed: number;
  experimentTypeId: string;
  experimentTypeName?: string;
  createdAt: number;
  updatedAt: number;
}

export interface ProtocolExperimentLink {
  id: Id<"protocol_experiment_links">;
  protocolId: Id<"protocols">;
  experimentId: string;
  createdAt: number;
  updatedAt: number;
}

export interface ProtocolExperimentResult {
  id: Id<"protocol_step_results">;
  protocolId: Id<"protocols">;
  experimentId: string;
  stepId: string;
  results: Record<string, any>;
  createdAt: number;
  updatedAt: number;
}

export interface CreateExperimentFromProtocolParams {
  title: string;
  description?: string;
  experimentTypeId: string;
  datePerformed?: number;
  parameters?: Record<string, any>;
  status?: string;
}

export class ProtocolIntegrationService {
  private convexClient: ConvexClient;
  
  constructor(convexClient: ConvexClient) {
    this.convexClient = convexClient;
  }
  
  /**
   * Get all experiments linked to a protocol
   * @param protocolId Protocol ID
   * @returns List of experiment summaries
   */
  async getLinkedExperiments(protocolId: Id<"protocols">): Promise<ExperimentSummary[]> {
    const links = await this.convexClient.query(api.protocolExperimentLinks.getLinksByProtocolId, { 
      protocolId 
    });
    
    if (!links || links.length === 0) {
      return [];
    }
    
    const experimentIds = links.map(link => link.experimentId);
    
    const experiments = await this.convexClient.query(api.experiments.enhancedExperiments.getExperimentsByIds, { 
      ids: experimentIds 
    });
    
    // Add experiment type names
    const experimentTypeIds = experiments.map(exp => exp.experimentTypeId).filter(Boolean);
    const experimentTypes = await this.convexClient.query(api.experiments.experimentTypes.getExperimentTypesByIds, {
      ids: experimentTypeIds
    });
    
    return experiments.map(exp => {
      const expType = experimentTypes.find(type => type._id === exp.experimentTypeId);
      
      return {
        id: exp._id,
        title: exp.title,
        description: exp.description,
        status: exp.status || 'planned',
        datePerformed: exp.datePerformed,
        experimentTypeId: exp.experimentTypeId,
        experimentTypeName: expType?.name,
        createdAt: exp._creationTime,
        updatedAt: exp._updatedAt || exp._creationTime,
      };
    });
  }
  
  /**
   * Get a protocol link by protocol ID and experiment ID
   * @param protocolId Protocol ID
   * @param experimentId Experiment ID
   * @returns Protocol-experiment link or null if not found
   */
  async getProtocolExperimentLink(
    protocolId: Id<"protocols">, 
    experimentId: string
  ): Promise<ProtocolExperimentLink | null> {
    const links = await this.convexClient.query(api.protocolExperimentLinks.getLinks, {
      protocolId,
      experimentId
    });
    
    if (!links || links.length === 0) {
      return null;
    }
    
    const link = links[0];
    
    return {
      id: link._id,
      protocolId: link.protocolId,
      experimentId: link.experimentId,
      createdAt: link._creationTime,
      updatedAt: link._updatedAt || link._creationTime,
    };
  }
  
  /**
   * Create a new experiment from a protocol
   * @param protocolId Protocol ID
   * @param params Experiment parameters
   * @returns Created experiment ID
   */
  async createExperimentFromProtocol(
    protocolId: Id<"protocols">,
    params: CreateExperimentFromProtocolParams
  ): Promise<string> {
    // 1. Create experiment
    const experimentId = await this.convexClient.mutation(api.experiments.enhancedExperiments.createExperiment, {
      title: params.title,
      description: params.description,
      experimentTypeId: params.experimentTypeId,
      datePerformed: params.datePerformed || Date.now(),
      parameters: params.parameters || {},
      status: params.status || "planned",
    });
    
    // 2. Create link between protocol and experiment
    await this.convexClient.mutation(api.protocolExperimentLinks.createLink, { 
      protocolId, 
      experimentId 
    });
    
    return experimentId;
  }
  
  /**
   * Delete a protocol-experiment link
   * @param linkId Link ID
   */
  async deleteLink(linkId: Id<"protocol_experiment_links">): Promise<void> {
    await this.convexClient.mutation(api.protocolExperimentLinks.deleteLink, { 
      linkId 
    });
  }
  
  /**
   * Get all results for an experiment related to a protocol
   * @param protocolId Protocol ID
   * @param experimentId Experiment ID
   * @returns List of step results
   */
  async getProtocolStepResults(
    protocolId: Id<"protocols">,
    experimentId: string
  ): Promise<ProtocolExperimentResult[]> {
    const results = await this.convexClient.query(api.protocolExperimentLinks.getAllResultsForExperiment, {
      experimentId
    });
    
    if (!results) {
      return [];
    }
    
    return results
      .filter(result => result.protocolId === protocolId)
      .map(result => ({
        id: result._id,
        protocolId: result.protocolId,
        experimentId: result.experimentId,
        stepId: result.stepId,
        results: result.results,
        createdAt: result._creationTime,
        updatedAt: result._updatedAt || result._creationTime,
      }));
  }
  
  /**
   * Record results for a protocol step in an experiment
   * @param protocolId Protocol ID
   * @param experimentId Experiment ID
   * @param stepId Step ID
   * @param results Results data
   * @returns Created/updated result ID
   */
  async recordStepResults(
    protocolId: Id<"protocols">,
    experimentId: string,
    stepId: string,
    results: Record<string, any>
  ): Promise<Id<"protocol_step_results">> {
    return await this.convexClient.mutation(api.protocolExperimentLinks.recordStepResults, {
      protocolId,
      experimentId,
      stepId,
      results
    });
  }
  
  /**
   * Get completion status for a protocol's steps in an experiment
   * @param protocolId Protocol ID
   * @param experimentId Experiment ID
   * @returns Completion information
   */
  async getProtocolCompletionStatus(
    protocolId: Id<"protocols">,
    experimentId: string
  ): Promise<{
    totalSteps: number;
    completedSteps: number;
    completionPercentage: number;
    stepStatus: Record<string, {
      completed: boolean;
      hasResults: boolean;
      lastUpdated?: number;
    }>;
  }> {
    const protocol = await this.convexClient.query(api.experiments.protocols.getProtocol, {
      protocolId
    });
    
    if (!protocol) {
      throw new Error(`Protocol with ID ${protocolId} not found`);
    }
    
    const results = await this.getProtocolStepResults(protocolId, experimentId);
    
    const stepStatus: Record<string, {
      completed: boolean;
      hasResults: boolean;
      lastUpdated?: number;
    }> = {};
    
    // Initialize step status for all steps
    protocol.steps.forEach(step => {
      stepStatus[step.id] = {
        completed: false,
        hasResults: false
      };
    });
    
    // Update step status based on results
    results.forEach(result => {
      if (stepStatus[result.stepId]) {
        stepStatus[result.stepId] = {
          completed: true,
          hasResults: true,
          lastUpdated: result.updatedAt
        };
      }
    });
    
    const totalSteps = protocol.steps.length;
    const completedSteps = Object.values(stepStatus).filter(status => status.completed).length;
    const completionPercentage = totalSteps > 0 ? Math.round((completedSteps / totalSteps) * 100) : 0;
    
    return {
      totalSteps,
      completedSteps,
      completionPercentage,
      stepStatus
    };
  }
}

export default ProtocolIntegrationService;