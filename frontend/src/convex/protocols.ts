/**
 * Protocol Template Functions for Convex
 * 
 * This file contains functions for managing protocol templates, including
 * creating, updating, exporting, and importing templates.
 */

import { v } from 'convex/values';
import { query, mutation } from './_generated/server';
import { Id } from './_generated/dataModel';

/**
 * Get a protocol template by ID
 */
export const getTemplate = query({
  args: {
    templateId: v.id('protocolTemplates')
  },
  handler: async (ctx, args) => {
    return await ctx.db.get(args.templateId);
  }
});

/**
 * List all protocol templates, optionally filtered by category
 */
export const listTemplates = query({
  args: {
    category: v.optional(v.string()),
    limit: v.optional(v.number()),
    excludeTemplateId: v.optional(v.id('protocolTemplates'))
  },
  handler: async (ctx, args) => {
    let q = ctx.db.query('protocolTemplates');
    
    if (args.category) {
      q = q.withIndex('by_category', q => q.eq('category', args.category));
    }
    
    // Exclude specific template if needed (for comparison view)
    if (args.excludeTemplateId) {
      q = q.filter(q => q.neq(q.field('_id'), args.excludeTemplateId));
    }
    
    // Sort by updatedAt in descending order (newest first)
    q = q.order('desc');
    
    // Apply limit if specified
    if (args.limit) {
      q = q.take(args.limit);
    }
    
    return await q.collect();
  }
});

/**
 * Create a new protocol template
 */
export const createTemplate = mutation({
  args: {
    name: v.string(),
    description: v.optional(v.string()),
    version: v.optional(v.string()),
    steps: v.array(
      v.object({
        name: v.string(),
        description: v.optional(v.string()),
        order: v.number(),
        duration: v.optional(v.number()),
        durationUnit: v.optional(v.string()),
        temperature: v.optional(v.number()),
        temperatureUnit: v.optional(v.string()),
        parameters: v.optional(v.map(v.string(), v.any())),
        equipment: v.optional(v.array(v.string())),
      })
    ),
    parameters: v.optional(v.map(v.string(), v.any())),
    category: v.optional(v.string()),
    tags: v.optional(v.array(v.string())),
    userId: v.optional(v.id('users'))
  },
  handler: async (ctx, args) => {
    // Default version if not provided
    const version = args.version || '1.0.0';
    
    // Create the template
    const templateId = await ctx.db.insert('protocolTemplates', {
      name: args.name,
      description: args.description,
      version: version,
      steps: args.steps,
      parameters: args.parameters || {},
      category: args.category,
      tags: args.tags || [],
      createdBy: args.userId,
      createdAt: Date.now(),
      updatedAt: Date.now()
    });
    
    return templateId;
  }
});

/**
 * Update an existing protocol template
 */
export const updateTemplate = mutation({
  args: {
    templateId: v.id('protocolTemplates'),
    name: v.optional(v.string()),
    description: v.optional(v.string()),
    steps: v.optional(v.array(
      v.object({
        name: v.string(),
        description: v.optional(v.string()),
        order: v.number(),
        duration: v.optional(v.number()),
        durationUnit: v.optional(v.string()),
        temperature: v.optional(v.number()),
        temperatureUnit: v.optional(v.string()),
        parameters: v.optional(v.map(v.string(), v.any())),
        equipment: v.optional(v.array(v.string())),
      })
    )),
    parameters: v.optional(v.map(v.string(), v.any())),
    category: v.optional(v.string()),
    tags: v.optional(v.array(v.string())),
    userId: v.optional(v.id('users'))
  },
  handler: async (ctx, args) => {
    // Get the existing template
    const template = await ctx.db.get(args.templateId);
    if (!template) {
      throw new Error('Template not found');
    }
    
    // Update the template
    await ctx.db.patch(args.templateId, {
      name: args.name !== undefined ? args.name : template.name,
      description: args.description !== undefined ? args.description : template.description,
      steps: args.steps !== undefined ? args.steps : template.steps,
      parameters: args.parameters !== undefined ? args.parameters : template.parameters,
      category: args.category !== undefined ? args.category : template.category,
      tags: args.tags !== undefined ? args.tags : template.tags,
      updatedAt: Date.now()
    });
    
    return args.templateId;
  }
});

/**
 * Create a new template version
 */
export const createTemplateVersion = mutation({
  args: {
    templateId: v.id('protocolTemplates'),
    name: v.optional(v.string()),
    description: v.optional(v.string()),
    version: v.string(),
    steps: v.optional(v.array(
      v.object({
        name: v.string(),
        description: v.optional(v.string()),
        order: v.number(),
        duration: v.optional(v.number()),
        durationUnit: v.optional(v.string()),
        temperature: v.optional(v.number()),
        temperatureUnit: v.optional(v.string()),
        parameters: v.optional(v.map(v.string(), v.any())),
        equipment: v.optional(v.array(v.string())),
      })
    )),
    parameters: v.optional(v.map(v.string(), v.any())),
    userId: v.optional(v.id('users'))
  },
  handler: async (ctx, args) => {
    // Get the parent template
    const parentTemplate = await ctx.db.get(args.templateId);
    if (!parentTemplate) {
      throw new Error('Parent template not found');
    }
    
    // Create the new version
    const newTemplateId = await ctx.db.insert('protocolTemplates', {
      name: args.name || parentTemplate.name,
      description: args.description !== undefined ? args.description : parentTemplate.description,
      version: args.version,
      parentVersionId: args.templateId,
      steps: args.steps || parentTemplate.steps,
      parameters: args.parameters || parentTemplate.parameters,
      category: parentTemplate.category,
      tags: parentTemplate.tags,
      createdBy: args.userId,
      createdAt: Date.now(),
      updatedAt: Date.now()
    });
    
    return newTemplateId;
  }
});

/**
 * Delete a protocol template
 */
export const deleteTemplate = mutation({
  args: {
    templateId: v.id('protocolTemplates')
  },
  handler: async (ctx, args) => {
    // Check if the template exists
    const template = await ctx.db.get(args.templateId);
    if (!template) {
      throw new Error('Template not found');
    }
    
    // Delete the template
    await ctx.db.delete(args.templateId);
    
    return args.templateId;
  }
});

/**
 * Export a protocol template to a portable format
 */
export const exportTemplate = mutation({
  args: {
    templateId: v.id('protocolTemplates')
  },
  handler: async (ctx, args) => {
    // Get the template
    const template = await ctx.db.get(args.templateId);
    if (!template) {
      throw new Error('Template not found');
    }
    
    // Format for export
    return {
      version: '1.0.0',
      exportDate: Date.now(),
      template: {
        id: args.templateId.toString(),
        name: template.name,
        description: template.description,
        version: template.version,
        parameters: template.parameters,
        steps: template.steps,
        category: template.category,
        tags: template.tags
      }
    };
  }
});

/**
 * Import a protocol template from a portable format
 */
export const importTemplate = mutation({
  args: {
    template: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      version: v.string(),
      parameters: v.optional(v.map(v.string(), v.any())),
      steps: v.array(
        v.object({
          name: v.string(),
          description: v.optional(v.string()),
          order: v.number(),
          duration: v.optional(v.number()),
          durationUnit: v.optional(v.string()),
          temperature: v.optional(v.number()),
          temperatureUnit: v.optional(v.string()),
          parameters: v.optional(v.map(v.string(), v.any())),
          equipment: v.optional(v.array(v.string())),
        })
      ),
      category: v.optional(v.string()),
      tags: v.optional(v.array(v.string()))
    }),
    userId: v.optional(v.id('users'))
  },
  handler: async (ctx, args) => {
    // Check if a template with the same name and version already exists
    const existingTemplates = await ctx.db
      .query('protocolTemplates')
      .withIndex('by_name_and_version', q => 
        q.eq('name', args.template.name)
         .eq('version', args.template.version)
      )
      .collect();
    
    if (existingTemplates.length > 0) {
      // Add a suffix to the version to make it unique
      args.template.version = `${args.template.version}-imported-${Date.now()}`;
    }
    
    // Create the imported template
    const templateId = await ctx.db.insert('protocolTemplates', {
      name: args.template.name,
      description: args.template.description,
      version: args.template.version,
      steps: args.template.steps,
      parameters: args.template.parameters || {},
      category: args.template.category,
      tags: args.template.tags || [],
      createdBy: args.userId,
      createdAt: Date.now(),
      updatedAt: Date.now()
    });
    
    return templateId;
  }
});