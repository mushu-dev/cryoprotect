/**
 * Service for integrating Protocol Designer with Experimental Data
 */

import { Protocol, ProtocolStep } from "../../experiments/services/experiment-service";
import { ProtocolService } from "./protocol-service";
import { ConvexClient } from "@convex/react";
import { api } from "../../../convex/_generated/api";
import { Id } from "../../../convex/_generated/dataModel";

export interface ExperimentProtocolLink {
  id: string;
  protocol_id: string;
  experiment_id: string;
  linked_at: number;
}

export interface ExperimentWithProtocol {
  id: string;
  title: string;
  experiment_type_id: string;
  description: string;
  protocol_id: string;
  protocol?: Protocol;
  date_performed?: string;
  temperature?: number;
  cooling_rate?: number;
  thawing_rate?: number;
  parameters?: Record<string, any>;
  version: number;
  created_at: string;
  updated_at: string;
  created_by: string;
}

export interface ProtocolWithExperiments {
  id: string;
  name: string;
  description: string;
  steps: ProtocolStep[];
  version: string;
  parent_id?: string;
  template_type: string;
  author?: string;
  parameters: Record<string, any>;
  experiments: any[];
}

export class ProtocolExperimentalService {
  private apiUrl: string;
  private convexClient: ConvexClient;

  constructor(baseUrl: string = '/api/v1', convexClient: ConvexClient) {
    this.apiUrl = baseUrl;
    this.convexClient = convexClient;
  }

  /**
   * Link an experiment with a protocol
   * @param experiment_id Experiment ID
   * @param protocol_id Protocol ID
   * @returns Link details
   */
  async linkExperimentToProtocol(experiment_id: string, protocol_id: string): Promise<ExperimentProtocolLink> {
    const response = await fetch(`${this.apiUrl}/experiment-protocol-links`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        experiment_id,
        protocol_id
      })
    });

    if (!response.ok) {
      throw new Error(`Failed to link experiment to protocol: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Get protocols linked to an experiment
   * @param experiment_id Experiment ID
   * @returns Array of links
   */
  async getProtocolsForExperiment(experiment_id: string): Promise<ExperimentProtocolLink[]> {
    const response = await fetch(`${this.apiUrl}/experiment-protocol-links?experiment_id=${experiment_id}`);

    if (!response.ok) {
      throw new Error(`Failed to get protocols for experiment: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Get experiments linked to a protocol
   * @param protocol_id Protocol ID
   * @returns Array of links
   */
  async getExperimentsForProtocol(protocol_id: string): Promise<ExperimentProtocolLink[]> {
    const response = await fetch(`${this.apiUrl}/experiment-protocol-links?protocol_id=${protocol_id}`);

    if (!response.ok) {
      throw new Error(`Failed to get experiments for protocol: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Delete a link between an experiment and protocol
   * @param link_id Link ID
   */
  async deleteLink(link_id: string): Promise<void> {
    const response = await fetch(`${this.apiUrl}/experiment-protocol-links/${link_id}`, {
      method: 'DELETE'
    });

    if (!response.ok) {
      throw new Error(`Failed to delete link: ${response.statusText}`);
    }
  }

  /**
   * Get an experiment with its linked protocol
   * @param experiment_id Experiment ID
   * @returns Experiment with protocol
   */
  async getExperimentWithProtocol(experiment_id: string): Promise<ExperimentWithProtocol> {
    const response = await fetch(`${this.apiUrl}/experiments/${experiment_id}/with-protocol`);

    if (!response.ok) {
      throw new Error(`Failed to get experiment with protocol: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Get a protocol with its linked experiments
   * @param protocol_id Protocol ID
   * @returns Protocol with experiments
   */
  async getProtocolWithExperiments(protocol_id: string): Promise<ProtocolWithExperiments> {
    const response = await fetch(`${this.apiUrl}/protocols/${protocol_id}/with-experiments`);

    if (!response.ok) {
      throw new Error(`Failed to get protocol with experiments: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Link a Convex protocol with an experiment using Convex API directly
   * This is an alternative implementation using Convex mutation directly
   * @param protocol_id Protocol ID (Convex ID)
   * @param experiment_id Experiment ID
   * @returns Link details
   */
  async linkProtocolToExperimentViaConvex(
    protocol_id: Id<"protocols">, 
    experiment_id: string
  ): Promise<Id<"protocol_experiment_links">> {
    return await this.convexClient.mutation(api.protocolExperimentLinks.createLink, { 
      protocolId: protocol_id,
      experimentId: experiment_id
    });
  }

  /**
   * Get results from experiment for a specific protocol step
   * @param protocol_id Protocol ID
   * @param step_id Step ID
   * @param experiment_id Experiment ID
   * @returns Results for the step
   */
  async getResultsForProtocolStep(
    protocol_id: string, 
    step_id: string, 
    experiment_id: string
  ): Promise<any> {
    const response = await fetch(
      `${this.apiUrl}/protocol-step-results?protocol_id=${protocol_id}&step_id=${step_id}&experiment_id=${experiment_id}`
    );

    if (!response.ok) {
      throw new Error(`Failed to get results for protocol step: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Record results for a protocol step in an experiment
   * @param protocol_id Protocol ID
   * @param step_id Step ID
   * @param experiment_id Experiment ID
   * @param results Results data
   * @returns Created result record
   */
  async recordStepResults(
    protocol_id: string,
    step_id: string,
    experiment_id: string,
    results: any
  ): Promise<any> {
    const response = await fetch(`${this.apiUrl}/protocol-step-results`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        protocol_id,
        step_id,
        experiment_id,
        results
      })
    });

    if (!response.ok) {
      throw new Error(`Failed to record step results: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data;
  }

  /**
   * Create an experiment from a protocol template
   * @param protocol_id Protocol ID
   * @param experiment_data Experiment data
   * @returns Created experiment ID
   */
  async createExperimentFromProtocol(
    protocol_id: string,
    experiment_data: {
      title: string;
      experiment_type_id: string;
      description?: string;
      date_performed?: string;
      parameters?: Record<string, any>;
    }
  ): Promise<string> {
    const response = await fetch(`${this.apiUrl}/experiments/from-protocol`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        protocol_id,
        experiment: experiment_data
      })
    });

    if (!response.ok) {
      throw new Error(`Failed to create experiment from protocol: ${response.statusText}`);
    }

    const result = await response.json();
    return result.data.id;
  }
}