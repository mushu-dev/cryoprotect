/**
 * CRUD operations for mixture components
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  MixtureComponent, 
  CreateMixtureComponentInput, 
  UpdateMixtureComponentInput,
  MixtureComponentFilter,
  MixtureComponentQueryOptions,
  MixtureComponentWithDetails 
} from "./types";
import { 
  validateCreateMixtureComponentInput, 
  validateUpdateMixtureComponentInput,
  validateMixtureComponentExists,
  validateMixtureExists,
  validateMixtureAccess
} from "./validation";
import { calculateTotalConcentration, validateConsistentUnits } from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new mixture component
 */
export const createMixtureComponent = mutation({
  args: {
    component: v.object({
      mixtureId: v.id("mixtures"),
      moleculeId: v.id("molecules"),
      concentration: v.number(),
      units: v.string(),
      role: v.optional(v.string()),
      notes: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateMixtureComponentInput(args.component);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get the mixture to check access
    const mixture = await ctx.db.get(args.component.mixtureId);
    if (!mixture) {
      throw new Error(`Mixture with ID ${args.component.mixtureId} not found`);
    }
    
    // Check access
    validateMixtureAccess(
      args.component.mixtureId, 
      userId, 
      mixture.public, 
      mixture.createdBy
    );
    
    // Check if the molecule exists
    const molecule = await ctx.db.get(args.component.moleculeId);
    if (!molecule) {
      throw new Error(`Molecule with ID ${args.component.moleculeId} not found`);
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const componentData = {
      ...args.component,
      createdAt: now,
      updatedAt: now
    };
    
    // Get existing components to check for consistency
    const existingComponents = await ctx.db
      .query("mixtureComponents")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.component.mixtureId))
      .collect();
    
    // Check for unit consistency if there are existing components
    if (existingComponents.length > 0) {
      const allComponents = [...existingComponents, componentData];
      if (!validateConsistentUnits(allComponents)) {
        throw new Error(
          "Inconsistent units: All components must use the same units"
        );
      }
    }
    
    // Insert component
    const componentId = await ctx.db.insert("mixtureComponents", componentData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "mixtureComponents",
      documentId: componentId,
      operation: "create",
      userId,
      newValue: componentData,
      timestamp: now,
    });
    
    // Update the mixture's updatedAt field
    await ctx.db.patch(args.component.mixtureId, { updatedAt: now });
    
    return componentId;
  }
});

/**
 * Get a mixture component by ID
 */
export const getMixtureComponent = query({
  args: {
    componentId: v.id("mixtureComponents"),
    includeMoleculeDetails: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate component exists
    validateMixtureComponentExists(args.componentId);
    
    // Get the component
    const component = await ctx.db.get(args.componentId);
    if (!component) {
      return null;
    }
    
    // Get the mixture to check access
    const mixture = await ctx.db.get(component.mixtureId);
    if (!mixture) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateMixtureAccess(
      component.mixtureId, 
      userId, 
      mixture.public, 
      mixture.createdBy
    );
    
    // Include molecule details if requested
    if (args.includeMoleculeDetails) {
      const molecule = await ctx.db.get(component.moleculeId);
      
      if (molecule) {
        return {
          ...component,
          molecule: {
            _id: molecule._id,
            name: molecule.name,
            formula: molecule.formula,
            pubchemCid: molecule.pubchemCid
          }
        };
      }
    }
    
    return component;
  }
});

/**
 * Update an existing mixture component
 */
export const updateMixtureComponent = mutation({
  args: {
    componentId: v.id("mixtureComponents"),
    update: v.object({
      concentration: v.optional(v.number()),
      units: v.optional(v.string()),
      role: v.optional(v.string()),
      notes: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateMixtureComponentInput(args.update);
    validateMixtureComponentExists(args.componentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing component
    const existingComponent = await ctx.db.get(args.componentId);
    if (!existingComponent) {
      throw new Error(`Mixture component with ID ${args.componentId} not found`);
    }
    
    // Get the mixture to check access
    const mixture = await ctx.db.get(existingComponent.mixtureId);
    if (!mixture) {
      throw new Error(`Mixture with ID ${existingComponent.mixtureId} not found`);
    }
    
    // Check access
    validateMixtureAccess(
      existingComponent.mixtureId, 
      userId, 
      mixture.public, 
      mixture.createdBy
    );
    
    // Check for unit consistency if changing units
    if (args.update.units && args.update.units !== existingComponent.units) {
      const existingComponents = await ctx.db
        .query("mixtureComponents")
        .withIndex("by_mixture", q => q.eq("mixtureId", existingComponent.mixtureId))
        .collect();
      
      // Create a temporary component with the new units
      const tempComponent = { 
        ...existingComponent, 
        units: args.update.units 
      };
      
      // Replace the current component with the temp one for checking
      const componentsToCheck = existingComponents
        .filter(c => !c._id.equals(args.componentId))
        .concat([tempComponent]);
      
      if (!validateConsistentUnits(componentsToCheck)) {
        throw new Error(
          "Inconsistent units: All components must use the same units"
        );
      }
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update component
    await ctx.db.patch(args.componentId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "mixtureComponents",
      documentId: args.componentId,
      operation: "update",
      userId,
      previousValue: existingComponent,
      newValue: { ...existingComponent, ...updateData },
      timestamp: now
    });
    
    // Update the mixture's updatedAt field
    await ctx.db.patch(existingComponent.mixtureId, { updatedAt: now });
    
    return args.componentId;
  }
});

/**
 * Delete a mixture component
 */
export const deleteMixtureComponent = mutation({
  args: {
    componentId: v.id("mixtureComponents")
  },
  handler: async (ctx, args) => {
    // Validate component exists
    validateMixtureComponentExists(args.componentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing component
    const existingComponent = await ctx.db.get(args.componentId);
    if (!existingComponent) {
      throw new Error(`Mixture component with ID ${args.componentId} not found`);
    }
    
    // Get the mixture to check access
    const mixture = await ctx.db.get(existingComponent.mixtureId);
    if (!mixture) {
      throw new Error(`Mixture with ID ${existingComponent.mixtureId} not found`);
    }
    
    // Check access
    validateMixtureAccess(
      existingComponent.mixtureId, 
      userId, 
      mixture.public, 
      mixture.createdBy
    );
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "mixtureComponents",
      documentId: args.componentId,
      operation: "delete",
      userId,
      previousValue: existingComponent,
      timestamp: now
    });
    
    // Delete the component
    await ctx.db.delete(args.componentId);
    
    // Update the mixture's updatedAt field
    await ctx.db.patch(existingComponent.mixtureId, { updatedAt: now });
    
    return true;
  }
});

/**
 * List components for a mixture
 */
export const listMixtureComponents = query({
  args: {
    mixtureId: v.id("mixtures"),
    includeMoleculeDetails: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate mixture exists
    validateMixtureExists(args.mixtureId);
    
    // Get the mixture to check access
    const mixture = await ctx.db.get(args.mixtureId);
    if (!mixture) {
      throw new Error(`Mixture with ID ${args.mixtureId} not found`);
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateMixtureAccess(
      args.mixtureId, 
      userId, 
      mixture.public, 
      mixture.createdBy
    );
    
    // Get components
    const components = await ctx.db
      .query("mixtureComponents")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.mixtureId))
      .collect();
    
    // Include molecule details if requested
    if (args.includeMoleculeDetails) {
      const componentsWithDetails: MixtureComponentWithDetails[] = [];
      
      for (const component of components) {
        const molecule = await ctx.db.get(component.moleculeId);
        
        if (molecule) {
          componentsWithDetails.push({
            ...component,
            molecule: {
              _id: molecule._id,
              name: molecule.name,
              formula: molecule.formula,
              pubchemCid: molecule.pubchemCid
            }
          });
        } else {
          componentsWithDetails.push({
            ...component,
            molecule: undefined
          });
        }
      }
      
      return componentsWithDetails;
    }
    
    return components;
  }
});

/**
 * Batch replace all components in a mixture
 */
export const replaceAllMixtureComponents = mutation({
  args: {
    mixtureId: v.id("mixtures"),
    components: v.array(
      v.object({
        moleculeId: v.id("molecules"),
        concentration: v.number(),
        units: v.string(),
        role: v.optional(v.string()),
        notes: v.optional(v.string())
      })
    )
  },
  handler: async (ctx, args) => {
    // Validate mixture exists
    validateMixtureExists(args.mixtureId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get the mixture to check access
    const mixture = await ctx.db.get(args.mixtureId);
    if (!mixture) {
      throw new Error(`Mixture with ID ${args.mixtureId} not found`);
    }
    
    // Check access - only creator can perform this operation
    if (!userId || !mixture.createdBy || !userId.equals(mixture.createdBy)) {
      throw new Error("Only the creator can replace all components");
    }
    
    // Validate input components
    for (const component of args.components) {
      validateCreateMixtureComponentInput({
        ...component,
        mixtureId: args.mixtureId
      });
      
      // Verify molecule exists
      const molecule = await ctx.db.get(component.moleculeId);
      if (!molecule) {
        throw new Error(`Molecule with ID ${component.moleculeId} not found`);
      }
    }
    
    // Check for unit consistency across all new components
    if (!validateConsistentUnits(
      args.components.map(c => ({ 
        ...c, 
        mixtureId: args.mixtureId,
        _id: undefined as any, 
        _creationTime: 0,
        createdAt: 0,
        updatedAt: 0
      }))
    )) {
      throw new Error(
        "Inconsistent units: All components must use the same units"
      );
    }
    
    // Get existing components
    const existingComponents = await ctx.db
      .query("mixtureComponents")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.mixtureId))
      .collect();
    
    const now = Date.now();
    
    // Delete all existing components with audit logs
    for (const component of existingComponents) {
      await ctx.db.insert("scientificDataAudit", {
        table: "mixtureComponents",
        documentId: component._id,
        operation: "delete",
        userId,
        previousValue: component,
        timestamp: now
      });
      
      await ctx.db.delete(component._id);
    }
    
    // Create all new components with audit logs
    const newComponentIds: Id<"mixtureComponents">[] = [];
    
    for (const component of args.components) {
      const componentData = {
        ...component,
        mixtureId: args.mixtureId,
        createdAt: now,
        updatedAt: now
      };
      
      const componentId = await ctx.db.insert("mixtureComponents", componentData);
      
      await ctx.db.insert("scientificDataAudit", {
        table: "mixtureComponents",
        documentId: componentId,
        operation: "create",
        userId,
        newValue: componentData,
        timestamp: now
      });
      
      newComponentIds.push(componentId);
    }
    
    // Update the mixture's updatedAt field
    await ctx.db.patch(args.mixtureId, { updatedAt: now });
    
    return newComponentIds;
  }
});