/**
 * Protocol Template Management System - Direct Convex Integration
 * This hook provides direct integration with Convex for protocol templates
 */

import { useQuery, useMutation } from "convex/react";
import { api } from "../../../../convex/_generated/api";
import { Id } from "../../../../convex/_generated/dataModel";
import { useErrorHandler } from "../../errors/use-error-handler";
import { useToast } from "../../ui/use-toast";

export type TemplateCategory = 
  | "Vitrification" 
  | "Slow Freezing" 
  | "Rapid Cooling" 
  | "Sample Preparation" 
  | "Thawing" 
  | "Analysis" 
  | "Custom";

export interface ProtocolTemplateFilter {
  name?: string;
  category?: TemplateCategory;
  public?: boolean;
  tags?: string[];
}

export interface ProtocolTemplateListOptions {
  limit?: number;
  cursor?: string;
  sortBy?: string;
  sortDirection?: "asc" | "desc";
}

export interface ProtocolTemplateCustomization {
  name: string;
  description?: string;
  steps?: {
    id: string;
    parameters?: Record<string, any>;
  }[];
  parameters?: Record<string, any>;
}

/**
 * Hook for managing protocol templates with direct Convex integration
 */
export function useProtocolTemplates(
  filter?: ProtocolTemplateFilter,
  options?: ProtocolTemplateListOptions
) {
  const { handleError } = useErrorHandler();
  const { toast } = useToast();

  // Query for listing protocol templates
  const templates = useQuery(
    api.experiments.protocols.listProtocols,
    filter ? { 
      filter: {
        ...filter,
        isTemplate: true 
      }, 
      options 
    } : { 
      filter: { isTemplate: true },
      options 
    }
  );

  // Get a single template
  const getTemplate = useQuery.load(api.experiments.protocols.getProtocol);

  // Create a new template
  const createTemplateMutation = useMutation(api.experiments.protocols.createProtocol);
  
  const createTemplate = async (template: any) => {
    try {
      const templateId = await createTemplateMutation({
        protocol: {
          ...template,
          isTemplate: true
        }
      });
      
      toast({
        title: "Template created",
        description: "Protocol template has been created successfully",
        variant: "success"
      });
      
      return templateId;
    } catch (error) {
      handleError(error, "Failed to create protocol template");
      throw error;
    }
  };

  // Update an existing template
  const updateTemplateMutation = useMutation(api.experiments.protocols.updateProtocol);
  
  const updateTemplate = async (templateId: Id<"protocols">, update: any) => {
    try {
      await updateTemplateMutation({
        protocolId: templateId,
        update
      });
      
      toast({
        title: "Template updated",
        description: "Protocol template has been updated successfully",
        variant: "success"
      });
      
      return templateId;
    } catch (error) {
      handleError(error, "Failed to update protocol template");
      throw error;
    }
  };

  // Delete a template
  const deleteTemplateMutation = useMutation(api.experiments.protocols.deleteProtocol);
  
  const deleteTemplate = async (templateId: Id<"protocols">) => {
    try {
      await deleteTemplateMutation({
        protocolId: templateId
      });
      
      toast({
        title: "Template deleted",
        description: "Protocol template has been deleted successfully",
        variant: "success"
      });
      
      return true;
    } catch (error) {
      handleError(error, "Failed to delete protocol template");
      throw error;
    }
  };

  // Create protocol from template
  const createProtocolFromTemplateMutation = useMutation(api.experiments.protocols.createProtocol);
  
  const createFromTemplate = async (
    templateId: Id<"protocols">, 
    customizations: ProtocolTemplateCustomization
  ) => {
    try {
      // First get the template
      const template = await getTemplate({ protocolId: templateId });
      
      if (!template) {
        throw new Error("Template not found");
      }
      
      // Apply customizations
      const customizedSteps = template.steps.map(step => {
        // Find if this step has customizations
        const stepCustomization = customizations.steps?.find(s => s.id === step._id.toString());
        
        if (stepCustomization) {
          return {
            ...step,
            parameters: {
              ...(step.parameters || {}),
              ...(stepCustomization.parameters || {})
            }
          };
        }
        
        return step;
      });
      
      // Create new protocol based on template
      const protocolId = await createProtocolFromTemplateMutation({
        protocol: {
          name: customizations.name,
          description: customizations.description || template.description,
          steps: customizedSteps,
          parameters: {
            ...(template.parameters || {}),
            ...(customizations.parameters || {})
          },
          parentId: templateId,
          isTemplate: false,
          category: template.category,
          public: false
        }
      });
      
      toast({
        title: "Protocol created",
        description: "Protocol has been created from template successfully",
        variant: "success"
      });
      
      return protocolId;
    } catch (error) {
      handleError(error, "Failed to create protocol from template");
      throw error;
    }
  };

  // Get template categories
  const getTemplateCategories = () => {
    return [
      "Vitrification",
      "Slow Freezing",
      "Rapid Cooling",
      "Sample Preparation",
      "Thawing",
      "Analysis",
      "Custom"
    ];
  };

  return {
    templates: templates || [],
    loading: templates === undefined,
    error: null,
    pagination: {
      hasMore: templates?.length === (options?.limit || 50)
    },
    getTemplate,
    createTemplate,
    updateTemplate,
    deleteTemplate,
    createFromTemplate,
    getTemplateCategories
  };
}

/**
 * Hook for accessing template comparison functionality
 */
export function useTemplateComparison() {
  const { handleError } = useErrorHandler();
  
  // Compare two templates
  const compareTemplates = async (
    templateId1: Id<"protocols">,
    templateId2: Id<"protocols">
  ) => {
    try {
      // Get both templates
      const template1 = await fetch(`/api/protocols/${templateId1}`).then(res => res.json());
      const template2 = await fetch(`/api/protocols/${templateId2}`).then(res => res.json());
      
      if (!template1 || !template2) {
        throw new Error("One or both templates not found");
      }
      
      // Compare the templates
      const comparison = {
        steps: {
          added: template2.steps.filter(step2 => !template1.steps.find(step1 => step1.name === step2.name)),
          removed: template1.steps.filter(step1 => !template2.steps.find(step2 => step2.name === step1.name)),
          modified: template1.steps
            .filter(step1 => {
              const step2 = template2.steps.find(s => s.name === step1.name);
              return step2 && JSON.stringify(step1) !== JSON.stringify(step2);
            })
            .map(step1 => {
              const step2 = template2.steps.find(s => s.name === step1.name);
              return {
                name: step1.name,
                changes: compareObjects(step1, step2)
              };
            })
        },
        parameters: compareObjects(template1.parameters || {}, template2.parameters || {})
      };
      
      return comparison;
    } catch (error) {
      handleError(error, "Failed to compare templates");
      throw error;
    }
  };
  
  // Helper function to compare objects
  const compareObjects = (obj1: any, obj2: any) => {
    const allKeys = [...new Set([...Object.keys(obj1), ...Object.keys(obj2)])];
    
    return {
      added: allKeys.filter(key => obj1[key] === undefined && obj2[key] !== undefined),
      removed: allKeys.filter(key => obj1[key] !== undefined && obj2[key] === undefined),
      modified: allKeys.filter(key => 
        obj1[key] !== undefined && 
        obj2[key] !== undefined && 
        JSON.stringify(obj1[key]) !== JSON.stringify(obj2[key])
      ).map(key => ({
        key,
        from: obj1[key],
        to: obj2[key]
      }))
    };
  };
  
  return {
    compareTemplates
  };
}

/**
 * Hook for managing template versions
 */
export function useTemplateVersions(templateName?: string) {
  const { handleError } = useErrorHandler();
  
  // Query for getting all versions of a template
  const versions = useQuery(
    api.experiments.protocols.getProtocolVersionHistory,
    templateName ? { name: templateName } : "skip"
  );
  
  return {
    versions: versions || [],
    loading: templateName !== undefined && versions === undefined,
    error: null
  };
}