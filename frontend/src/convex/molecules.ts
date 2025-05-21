import { query } from "../../../convex/_generated/server";
import { v } from "convex/values";

/**
 * Query to get all molecules from Convex
 */
export const getAllMolecules = query({
  args: {
    limit: v.optional(v.number())
  },
  async handler(ctx, args) {
    // Set a reasonable default limit
    const limit = args.limit || 100;
    
    // Query molecules from the database
    const molecules = await ctx.db
      .query("molecules")
      .order("desc")
      .take(limit);
    
    // For each molecule, fetch its properties
    const moleculesWithProperties = await Promise.all(
      molecules.map(async (molecule) => {
        // Get properties for this molecule
        const properties = await ctx.db
          .query("molecularProperties")
          .filter(q => q.eq("moleculeId", molecule._id))
          .collect();
        
        // Process properties to make them more usable
        const processedProperties = {};
        for (const property of properties) {
          // Get property type details
          const propertyType = await ctx.db.get(property.propertyTypeId);
          if (propertyType) {
            // Use property type name as the key
            processedProperties[propertyType.name] = {
              value: property.value,
              units: property.units || propertyType.units,
              displayName: propertyType.displayName
            };
          }
        }
        
        // Return molecule with properties
        return {
          id: molecule._id,
          name: molecule.name,
          formula: molecule.formula,
          smiles: molecule.canonicalSmiles,
          pubchemCid: molecule.pubchemCid,
          properties: processedProperties
        };
      })
    );
    
    return moleculesWithProperties;
  }
});

/**
 * Query to get a single molecule by ID from Convex
 */
export const getMoleculeById = query({
  args: {
    id: v.string()
  },
  async handler(ctx, args) {
    // Get the molecule by ID
    const molecule = await ctx.db.get(args.id);
    
    if (!molecule) {
      return null;
    }
    
    // Get properties for this molecule
    const properties = await ctx.db
      .query("molecularProperties")
      .filter(q => q.eq("moleculeId", molecule._id))
      .collect();
    
    // Process properties
    const processedProperties = {};
    for (const property of properties) {
      // Get property type details
      const propertyType = await ctx.db.get(property.propertyTypeId);
      if (propertyType) {
        // Use property type name as the key
        processedProperties[propertyType.name] = {
          value: property.value,
          units: property.units || propertyType.units,
          displayName: propertyType.displayName
        };
      }
    }
    
    // Return molecule with properties
    return {
      id: molecule._id,
      name: molecule.name,
      formula: molecule.formula,
      smiles: molecule.canonicalSmiles,
      pubchemCid: molecule.pubchemCid,
      properties: processedProperties
    };
  }
});