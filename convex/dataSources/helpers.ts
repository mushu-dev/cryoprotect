/**
 * Helper functions for data sources
 */

import { Id } from "../_generated/dataModel";
import { QueryCtx } from "../_generated/server";
import { DataSource, DataSourceWithUsage } from "./types";

/**
 * Format data source type for display
 */
export function formatDataSourceType(type: string): string {
  switch (type) {
    case "database":
      return "Database";
    case "publication":
      return "Publication";
    case "experiment":
      return "Experiment";
    case "calculation":
      return "Calculation";
    default:
      return type.charAt(0).toUpperCase() + type.slice(1);
  }
}

/**
 * Get statistics for a data source's usage in the system
 */
export async function getDataSourceUsageStats(
  ctx: QueryCtx,
  dataSourceId: Id<"dataSources">
): Promise<{ moleculeCount: number; propertyCount: number; experimentCount: number }> {
  // Count molecules using this data source
  const moleculeCount = await ctx.db
    .query("molecules")
    .filter(q => q.eq(q.field("dataSource"), dataSourceId))
    .count();
  
  // Count properties using this data source
  const propertyCount = await ctx.db
    .query("molecularProperties")
    .filter(q => q.eq(q.field("source"), dataSourceId))
    .count();
  
  // Count experiments referencing this data source (assuming such a link exists)
  // In a real implementation, we would have a specific field for this
  // For now, we'll always return 0
  const experimentCount = 0;
  
  return {
    moleculeCount,
    propertyCount,
    experimentCount
  };
}

/**
 * Expand a data source with its usage statistics
 */
export async function expandDataSourceWithUsage(
  ctx: QueryCtx,
  dataSource: DataSource
): Promise<DataSourceWithUsage> {
  const usageStats = await getDataSourceUsageStats(ctx, dataSource._id);
  
  return {
    ...dataSource,
    usageStats
  };
}

/**
 * Format publication citation from data source
 */
export function formatPublicationCitation(dataSource: DataSource): string {
  if (dataSource.type !== "publication") {
    return dataSource.name;
  }
  
  // In a real implementation, we would have more specific fields for publications
  // such as authors, journal, year, etc.
  // For now, we'll just return the name and URL if available
  if (dataSource.url) {
    return `${dataSource.name}. ${dataSource.url}`;
  }
  
  return dataSource.name;
}

/**
 * Check if a data source is in use
 */
export async function isDataSourceInUse(
  ctx: QueryCtx,
  dataSourceId: Id<"dataSources">
): Promise<boolean> {
  // Check if any molecules use this data source
  const moleculeCount = await ctx.db
    .query("molecules")
    .filter(q => q.eq(q.field("dataSource"), dataSourceId))
    .count();
  
  if (moleculeCount > 0) {
    return true;
  }
  
  // Check if any properties use this data source
  const propertyCount = await ctx.db
    .query("molecularProperties")
    .filter(q => q.eq(q.field("source"), dataSourceId))
    .count();
  
  if (propertyCount > 0) {
    return true;
  }
  
  // In a real implementation, we would check other tables as well
  
  return false;
}