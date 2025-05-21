/**
 * Helper functions for mixtures and mixture components
 */

import { Id } from "../_generated/dataModel";
import { DataModel } from "../_generated/dataModel";
import { ConvexError } from "convex/values";
import { QueryCtx } from "../_generated/server";
import { Mixture, MixtureComponent, MixtureWithComponents, MixtureComponentWithDetails } from "./types";

/**
 * Calculates the total concentration in a mixture to ensure valid proportions
 */
export function calculateTotalConcentration(components: MixtureComponent[]): number {
  if (!components || components.length === 0) {
    return 0;
  }
  
  // Group components by units
  const concentrationsByUnit: Record<string, number> = {};
  
  for (const component of components) {
    if (!concentrationsByUnit[component.units]) {
      concentrationsByUnit[component.units] = 0;
    }
    concentrationsByUnit[component.units] += component.concentration;
  }
  
  // For simplicity, we'll return the sum of the first unit type found
  // In a real implementation, we would do unit conversion
  const firstUnit = Object.keys(concentrationsByUnit)[0];
  return concentrationsByUnit[firstUnit] || 0;
}

/**
 * Validates component units are consistent
 * Returns true if all components use the same units, false otherwise
 */
export function validateConsistentUnits(components: MixtureComponent[]): boolean {
  if (!components || components.length <= 1) {
    return true;
  }
  
  const unitSet = new Set(components.map(c => c.units));
  return unitSet.size === 1;
}

/**
 * Normalizes mixture component concentrations to sum to 100%
 * Only applicable for percentages or similar unit types
 */
export function normalizeConcentrations(
  components: MixtureComponent[]
): MixtureComponent[] {
  if (!components || components.length === 0) {
    return components;
  }
  
  // Check if normalization is applicable (only for percentage-like units)
  const percentageUnits = ["percent", "percentage", "%", "w/w", "v/v", "w/v"];
  const areAllPercentage = components.every(c => 
    percentageUnits.includes(c.units.toLowerCase())
  );
  
  if (!areAllPercentage) {
    return components;
  }
  
  // Calculate total
  const total = components.reduce((sum, c) => sum + c.concentration, 0);
  
  // If total is already 100 (or very close), return as is
  if (Math.abs(total - 100) < 0.001) {
    return components;
  }
  
  // Normalize to 100%
  return components.map(c => ({
    ...c,
    concentration: (c.concentration / total) * 100
  }));
}

/**
 * Expand a mixture with its components
 */
export async function expandMixtureWithComponents(
  ctx: QueryCtx,
  mixture: Mixture
): Promise<MixtureWithComponents> {
  const components = await ctx.db
    .query("mixtureComponents")
    .withIndex("by_mixture", q => q.eq("mixtureId", mixture._id))
    .collect();
  
  const componentWithDetails: MixtureComponentWithDetails[] = [];
  
  for (const component of components) {
    const molecule = await ctx.db.get(component.moleculeId);
    componentWithDetails.push({
      ...component,
      molecule: molecule ? {
        _id: molecule._id,
        name: molecule.name,
        formula: molecule.formula,
        pubchemCid: molecule.pubchemCid
      } : undefined
    });
  }
  
  return {
    ...mixture,
    components: componentWithDetails
  };
}

/**
 * Format concentration with appropriate units
 */
export function formatConcentration(concentration: number, units: string): string {
  // Handle percentage units
  if (units === "%" || units.toLowerCase().includes("percent")) {
    return `${concentration.toFixed(2)}${units}`;
  }
  
  // Handle molar concentrations
  if (units.toLowerCase().includes("mol") || 
      units.toLowerCase().includes("m") ||
      units.toLowerCase() === "mm" ||
      units.toLowerCase() === "μm") {
    return `${concentration.toFixed(4)} ${units}`;
  }
  
  // Handle weight/volume
  if (units.toLowerCase().includes("w/v") || 
      units.toLowerCase().includes("mg/ml") ||
      units.toLowerCase().includes("g/l")) {
    return `${concentration.toFixed(2)} ${units}`;
  }
  
  // Default formatting
  return `${concentration} ${units}`;
}