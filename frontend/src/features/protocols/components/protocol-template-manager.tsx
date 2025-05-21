/**
 * Protocol Template Manager Component
 * 
 * A comprehensive interface for managing protocol templates, including
 * browsing, creating, editing, and importing/exporting templates.
 */

import React, { useState } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';
import { Id } from '../../../convex/_generated/dataModel';
import { TemplateExportImport } from './template-export-import';
import { useToast } from '../../ui/use-toast';

/**
 * Protocol Template Manager Component
 */
export function ProtocolTemplateManager() {
  // State for selected template
  const [selectedTemplateId, setSelectedTemplateId] = useState<Id<"protocolTemplates"> | null>(null);
  
  // State for view mode
  const [viewMode, setViewMode] = useState<'list' | 'detail' | 'export-import'>('list');
  
  // State for search and filter
  const [searchTerm, setSearchTerm] = useState('');
  const [selectedCategory, setSelectedCategory] = useState<string | null>(null);
  
  // Toast for notifications
  const { toast } = useToast();
  
  // Get all templates
  const templates = useQuery(api.protocols.templates.listTemplates, {});
  
  // Get categories from templates
  const categories = templates
    ? Array.from(new Set(templates.filter(t => t.category).map(t => t.category)))
    : [];
  
  // Filter templates based on search and category
  const filteredTemplates = templates
    ? templates.filter(template => {
        const matchesSearch = searchTerm
          ? template.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
            (template.description && template.description.toLowerCase().includes(searchTerm.toLowerCase()))
          : true;
        
        const matchesCategory = selectedCategory
          ? template.category === selectedCategory
          : true;
        
        return matchesSearch && matchesCategory;
      })
    : [];
  
  // Get selected template details
  const selectedTemplate = useQuery(
    api.protocols.templates.getTemplate,
    selectedTemplateId ? { templateId: selectedTemplateId } : "skip"
  );
  
  // Delete template mutation
  const deleteTemplate = useMutation(api.protocols.templates.deleteTemplate);
  
  // Handle template selection
  const handleSelectTemplate = (templateId: Id<"protocolTemplates">) => {
    setSelectedTemplateId(templateId);
    setViewMode('detail');
  };
  
  // Handle template deletion
  const handleDeleteTemplate = async (templateId: Id<"protocolTemplates">) => {
    if (window.confirm('Are you sure you want to delete this template? This action cannot be undone.')) {
      try {
        await deleteTemplate({ templateId });
        toast({
          title: 'Template Deleted',
          description: 'The protocol template has been deleted successfully.',
        });
        
        if (selectedTemplateId === templateId) {
          setSelectedTemplateId(null);
          setViewMode('list');
        }
      } catch (error) {
        toast({
          title: 'Delete Failed',
          description: 'Failed to delete the template. Please try again.',
          variant: 'destructive'
        });
      }
    }
  };
  
  // Handle import completion
  const handleImportComplete = (importedTemplateId: Id<"protocolTemplates">) => {
    setSelectedTemplateId(importedTemplateId);
    setViewMode('detail');
    toast({
      title: 'Template Imported',
      description: 'Protocol template has been imported successfully.',
    });
  };
  
  return (
    <div className="container mx-auto py-6">
      <div className="flex justify-between items-center mb-6">
        <h1 className="text-2xl font-bold">Protocol Templates</h1>
        
        <div className="flex space-x-4">
          <button
            onClick={() => setViewMode('export-import')}
            className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
          >
            Import / Export
          </button>
          <button
            onClick={() => {
              // Handle creating a new template
              // This would typically navigate to a template creation page
              alert('Navigate to template creation page');
            }}
            className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 transition-colors"
          >
            Create Template
          </button>
        </div>
      </div>
      
      {viewMode === 'list' && (
        <>
          {/* Search and filter */}
          <div className="mb-6 flex flex-col md:flex-row md:items-center space-y-3 md:space-y-0 md:space-x-4">
            <div className="flex-1">
              <input
                type="text"
                value={searchTerm}
                onChange={(e) => setSearchTerm(e.target.value)}
                placeholder="Search templates..."
                className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
              />
            </div>
            
            <div>
              <select
                value={selectedCategory || ''}
                onChange={(e) => setSelectedCategory(e.target.value || null)}
                className="px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
              >
                <option value="">All Categories</option>
                {categories.map((category, index) => (
                  <option key={index} value={category}>{category}</option>
                ))}
              </select>
            </div>
          </div>
          
          {/* Templates list */}
          {templates === undefined ? (
            <div className="text-center py-8">Loading templates...</div>
          ) : templates.length === 0 ? (
            <div className="text-center py-8">
              <p className="text-gray-500">No templates found.</p>
              <p className="text-gray-500 mt-2">Create a new template or import one to get started.</p>
            </div>
          ) : filteredTemplates.length === 0 ? (
            <div className="text-center py-8">
              <p className="text-gray-500">No templates match your search criteria.</p>
            </div>
          ) : (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {filteredTemplates.map((template) => (
                <div 
                  key={template._id.toString()}
                  className="border border-gray-300 rounded-lg overflow-hidden hover:shadow-md transition-shadow"
                >
                  <div className="px-5 py-4">
                    <div className="flex justify-between items-start">
                      <h3 
                        className="text-lg font-semibold mb-2 cursor-pointer text-blue-600 hover:text-blue-800"
                        onClick={() => handleSelectTemplate(template._id)}
                      >
                        {template.name}
                      </h3>
                      <span className="bg-gray-200 text-gray-700 px-2 py-1 rounded-full text-xs">
                        v{template.version}
                      </span>
                    </div>
                    
                    {template.description && (
                      <p className="text-gray-600 text-sm mb-3 line-clamp-2">
                        {template.description}
                      </p>
                    )}
                    
                    <div className="flex justify-between items-center text-xs text-gray-500">
                      <span>{template.steps.length} steps</span>
                      <span>
                        {new Date(template.updatedAt).toLocaleDateString()}
                      </span>
                    </div>
                    
                    {template.category && (
                      <span className="inline-block mt-2 bg-blue-100 text-blue-800 px-2 py-1 rounded-full text-xs">
                        {template.category}
                      </span>
                    )}
                    
                    {template.tags && template.tags.length > 0 && (
                      <div className="mt-3 flex flex-wrap gap-1">
                        {template.tags.slice(0, 3).map((tag, index) => (
                          <span 
                            key={index}
                            className="px-2 py-0.5 bg-gray-100 text-gray-700 rounded-full text-xs"
                          >
                            {tag}
                          </span>
                        ))}
                        {template.tags.length > 3 && (
                          <span className="px-2 py-0.5 bg-gray-100 text-gray-700 rounded-full text-xs">
                            +{template.tags.length - 3} more
                          </span>
                        )}
                      </div>
                    )}
                  </div>
                  
                  <div className="bg-gray-50 px-5 py-3 flex justify-end space-x-2">
                    <button
                      onClick={() => handleSelectTemplate(template._id)}
                      className="px-3 py-1 text-sm bg-transparent text-blue-600 hover:text-blue-800"
                    >
                      View
                    </button>
                    <button
                      onClick={() => {
                        // Handle editing the template
                        alert(`Edit template ${template._id}`);
                      }}
                      className="px-3 py-1 text-sm bg-transparent text-gray-600 hover:text-gray-800"
                    >
                      Edit
                    </button>
                    <button
                      onClick={() => handleDeleteTemplate(template._id)}
                      className="px-3 py-1 text-sm bg-transparent text-red-600 hover:text-red-800"
                    >
                      Delete
                    </button>
                  </div>
                </div>
              ))}
            </div>
          )}
        </>
      )}
      
      {viewMode === 'detail' && selectedTemplate && (
        <div className="bg-white border border-gray-300 rounded-lg overflow-hidden">
          <div className="px-6 py-5 border-b border-gray-300 flex justify-between items-center">
            <button
              onClick={() => {
                setSelectedTemplateId(null);
                setViewMode('list');
              }}
              className="text-blue-600 hover:text-blue-800"
            >
              ← Back to templates
            </button>
            
            <div className="flex space-x-3">
              <button
                onClick={() => setViewMode('export-import')}
                className="px-3 py-1 text-sm bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Export
              </button>
              <button
                onClick={() => {
                  // Handle editing the template
                  alert(`Edit template ${selectedTemplate._id}`);
                }}
                className="px-3 py-1 text-sm bg-gray-200 text-gray-800 rounded-md hover:bg-gray-300"
              >
                Edit
              </button>
            </div>
          </div>
          
          <div className="px-6 py-5">
            <div className="flex justify-between items-start mb-4">
              <h2 className="text-2xl font-bold">{selectedTemplate.name}</h2>
              <span className="bg-gray-200 text-gray-700 px-3 py-1 rounded-full text-sm">
                v{selectedTemplate.version}
              </span>
            </div>
            
            {selectedTemplate.description && (
              <p className="text-gray-700 mb-6">
                {selectedTemplate.description}
              </p>
            )}
            
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
              <div className="bg-gray-50 p-4 rounded-md">
                <h3 className="text-sm font-medium text-gray-500 mb-1">Category</h3>
                <p>{selectedTemplate.category || 'Uncategorized'}</p>
              </div>
              
              <div className="bg-gray-50 p-4 rounded-md">
                <h3 className="text-sm font-medium text-gray-500 mb-1">Steps</h3>
                <p>{selectedTemplate.steps.length} steps</p>
              </div>
              
              <div className="bg-gray-50 p-4 rounded-md">
                <h3 className="text-sm font-medium text-gray-500 mb-1">Last Updated</h3>
                <p>{new Date(selectedTemplate.updatedAt).toLocaleString()}</p>
              </div>
            </div>
            
            {selectedTemplate.tags && selectedTemplate.tags.length > 0 && (
              <div className="mb-6">
                <h3 className="text-sm font-medium text-gray-500 mb-2">Tags</h3>
                <div className="flex flex-wrap gap-2">
                  {selectedTemplate.tags.map((tag, index) => (
                    <span 
                      key={index}
                      className="px-3 py-1 bg-blue-100 text-blue-800 rounded-full text-sm"
                    >
                      {tag}
                    </span>
                  ))}
                </div>
              </div>
            )}
            
            <div className="mb-6">
              <h3 className="text-lg font-medium mb-4">Protocol Steps</h3>
              
              <div className="space-y-4">
                {selectedTemplate.steps
                  .sort((a, b) => a.order - b.order)
                  .map((step, index) => (
                    <div 
                      key={index}
                      className="border border-gray-200 rounded-md overflow-hidden"
                    >
                      <div className="bg-gray-50 px-4 py-3 flex justify-between items-center">
                        <div className="flex items-center">
                          <span className="inline-flex items-center justify-center bg-blue-600 text-white rounded-full w-6 h-6 text-sm font-medium mr-3">
                            {step.order}
                          </span>
                          <h4 className="font-medium">{step.name}</h4>
                        </div>
                        
                        {step.duration && step.durationUnit && (
                          <span className="text-sm text-gray-600">
                            {step.duration} {step.durationUnit}
                          </span>
                        )}
                      </div>
                      
                      {step.description && (
                        <div className="px-4 py-3">
                          <p className="text-gray-700 text-sm">{step.description}</p>
                        </div>
                      )}
                      
                      <div className="px-4 py-3 border-t border-gray-200 grid grid-cols-1 md:grid-cols-2 gap-4">
                        {step.temperature && (
                          <div>
                            <span className="text-xs text-gray-500">Temperature</span>
                            <p className="text-sm">
                              {step.temperature} {step.temperatureUnit || '°C'}
                            </p>
                          </div>
                        )}
                        
                        {step.equipment && step.equipment.length > 0 && (
                          <div>
                            <span className="text-xs text-gray-500">Equipment</span>
                            <p className="text-sm">
                              {step.equipment.join(', ')}
                            </p>
                          </div>
                        )}
                      </div>
                    </div>
                  ))}
              </div>
            </div>
            
            {/* Parameters section */}
            {selectedTemplate.parameters && Object.keys(selectedTemplate.parameters).length > 0 && (
              <div className="mb-6">
                <h3 className="text-lg font-medium mb-4">Parameters</h3>
                
                <div className="border border-gray-200 rounded-md overflow-hidden">
                  <table className="min-w-full divide-y divide-gray-200">
                    <thead className="bg-gray-50">
                      <tr>
                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                          Parameter
                        </th>
                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                          Value
                        </th>
                      </tr>
                    </thead>
                    <tbody className="bg-white divide-y divide-gray-200">
                      {Object.entries(selectedTemplate.parameters).map(([key, value], index) => (
                        <tr key={index}>
                          <td className="px-4 py-3 text-sm text-gray-900">{key}</td>
                          <td className="px-4 py-3 text-sm text-gray-500">
                            {typeof value === 'object' 
                              ? JSON.stringify(value) 
                              : String(value)}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}
            
            <div className="flex justify-between pt-4 border-t border-gray-200">
              <button
                onClick={() => {
                  // Clone the template as a starting point for a new one
                  alert(`Clone template ${selectedTemplate._id}`);
                }}
                className="px-4 py-2 bg-gray-200 text-gray-800 rounded-md hover:bg-gray-300 transition-colors"
              >
                Clone Template
              </button>
              
              <button
                onClick={() => {
                  // Create a protocol from this template
                  alert(`Create protocol from template ${selectedTemplate._id}`);
                }}
                className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 transition-colors"
              >
                Create Protocol from Template
              </button>
            </div>
          </div>
        </div>
      )}
      
      {viewMode === 'export-import' && (
        <div className="bg-white border border-gray-300 rounded-lg overflow-hidden">
          <div className="px-6 py-5 border-b border-gray-300">
            <button
              onClick={() => {
                if (selectedTemplateId) {
                  setViewMode('detail');
                } else {
                  setViewMode('list');
                }
              }}
              className="text-blue-600 hover:text-blue-800"
            >
              ← {selectedTemplateId ? 'Back to template details' : 'Back to templates'}
            </button>
          </div>
          
          <div className="p-6">
            <h2 className="text-2xl font-bold mb-6">Import/Export Protocol Templates</h2>
            
            <TemplateExportImport
              templateId={selectedTemplateId || undefined}
              onImportComplete={handleImportComplete}
              onExportComplete={() => {
                toast({
                  title: 'Template Exported',
                  description: 'Protocol template has been exported successfully.',
                });
              }}
            />
          </div>
        </div>
      )}
    </div>
  );
}