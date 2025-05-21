/**
 * Type definitions for data sources
 */

import { Id } from "../_generated/dataModel";

/**
 * Data source document in the database
 */
export interface DataSource {
  _id: Id<"dataSources">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  url?: string;
  type: "database" | "publication" | "experiment" | "calculation";
  version?: string;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new data source
 */
export interface CreateDataSourceInput {
  name: string;
  description?: string;
  url?: string;
  type: "database" | "publication" | "experiment" | "calculation";
  version?: string;
}

/**
 * Input for updating a data source
 */
export interface UpdateDataSourceInput {
  name?: string;
  description?: string;
  url?: string;
  type?: "database" | "publication" | "experiment" | "calculation";
  version?: string;
}

/**
 * Data source with usage statistics
 */
export interface DataSourceWithUsage extends DataSource {
  usageStats: {
    moleculeCount: number;
    propertyCount: number;
    experimentCount: number;
  };
}

/**
 * Filter conditions for querying data sources
 */
export interface DataSourceFilter {
  name?: string;
  type?: "database" | "publication" | "experiment" | "calculation";
  version?: string;
}

/**
 * Options for data source queries
 */
export interface DataSourceQueryOptions {
  limit?: number;
  cursor?: string;
  includeUsageStats?: boolean;
  sortBy?: "name" | "type" | "updatedAt";
  sortDirection?: "asc" | "desc";
}