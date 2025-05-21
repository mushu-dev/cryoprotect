/**
 * Protocol Template Export/Import Component
 * 
 * This component allows users to export protocol templates to JSON files
 * and import them back into the system. It supports versioning, validation,
 * and compatibility checking.
 */

import React, { useState, useRef } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';
import { Id } from '../../../convex/_generated/dataModel';
import { useToast } from '../../ui/use-toast';
import { useErrorHandler } from '../../errors/use-error-handler';

// Types
export interface ProtocolTemplateExportData {
  version: string;
  exportDate: number;
  template: {
    id: string;
    name: string;
    description?: string;
    version: string;
    parameters?: Record<string, any>;
    steps: Array<{
      name: string;
      description?: string;
      order: number;
      duration?: number;
      durationUnit?: string;
      temperature?: number;
      temperatureUnit?: string;
      parameters?: Record<string, any>;
      equipment?: string[];
    }>;
    category?: string;
    tags?: string[];
  };
}

interface TemplateExportImportProps {
  templateId?: Id<"protocolTemplates">;
  onImportComplete?: (importedTemplateId: Id<"protocolTemplates">) => void;
  onExportComplete?: () => void;
  showExport?: boolean;
  showImport?: boolean;
  className?: string;
}

/**
 * Protocol Template Export/Import Component
 */
export function TemplateExportImport({
  templateId,
  onImportComplete,
  onExportComplete,
  showExport = true,
  showImport = true,
  className = ''
}: TemplateExportImportProps) {
  // File input ref
  const fileInputRef = useRef<HTMLInputElement>(null);
  
  // Import states
  const [isImporting, setIsImporting] = useState(false);
  const [importFile, setImportFile] = useState<File | null>(null);
  const [importError, setImportError] = useState<string | null>(null);
  const [importData, setImportData] = useState<ProtocolTemplateExportData | null>(null);
  const [showImportPreview, setShowImportPreview] = useState(false);
  
  // Export states
  const [isExporting, setIsExporting] = useState(false);
  
  // Toast notifications
  const { toast } = useToast();
  const { handleError } = useErrorHandler();
  
  // Get the template data if templateId is provided
  const template = useQuery(
    api.protocols.templates.getTemplate,
    templateId ? { templateId } : "skip"
  );
  
  // Mutations
  const importTemplateMutation = useMutation(api.protocols.templates.importTemplate);
  const exportTemplateMutation = useMutation(api.protocols.templates.exportTemplate);
  
  // Handle file selection for import
  const handleFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    setImportError(null);
    setImportData(null);
    setShowImportPreview(false);
    
    const files = e.target.files;
    if (!files || files.length === 0) {
      return;
    }
    
    const file = files[0];
    setImportFile(file);
    
    // Read the file contents
    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        if (event.target && typeof event.target.result === 'string') {
          const data = JSON.parse(event.target.result) as ProtocolTemplateExportData;
          
          // Validate the imported data
          if (!validateImportData(data)) {
            setImportError('Invalid template format. Please check the file and try again.');
            return;
          }
          
          setImportData(data);
          setShowImportPreview(true);
        }
      } catch (error) {
        setImportError('Failed to parse template file. Please check the file format and try again.');
        console.error('Import error:', error);
      }
    };
    
    reader.onerror = () => {
      setImportError('Failed to read the file. Please try again.');
    };
    
    reader.readAsText(file);
  };
  
  // Validate imported data format
  const validateImportData = (data: any): data is ProtocolTemplateExportData => {
    if (!data.version || !data.exportDate || !data.template) {
      return false;
    }
    
    const template = data.template;
    if (!template.name || !template.version || !Array.isArray(template.steps)) {
      return false;
    }
    
    // Validate each step has required fields
    for (const step of template.steps) {
      if (!step.name || typeof step.order !== 'number') {
        return false;
      }
    }
    
    return true;
  };
  
  // Handle template import
  const handleImport = async () => {
    if (!importData) {
      setImportError('No template data to import.');
      return;
    }
    
    setIsImporting(true);
    
    try {
      const templateData = {
        name: importData.template.name,
        description: importData.template.description,
        version: importData.template.version,
        parameters: importData.template.parameters || {},
        steps: importData.template.steps.map(step => ({
          name: step.name,
          description: step.description,
          order: step.order,
          duration: step.duration,
          durationUnit: step.durationUnit,
          temperature: step.temperature,
          temperatureUnit: step.temperatureUnit,
          parameters: step.parameters || {},
          equipment: step.equipment || []
        })),
        category: importData.template.category,
        tags: importData.template.tags || []
      };
      
      const importedTemplateId = await importTemplateMutation({ template: templateData });
      
      toast({
        title: 'Template Imported',
        description: `Successfully imported template "${importData.template.name}".`,
      });
      
      // Reset states
      setImportFile(null);
      setImportData(null);
      setShowImportPreview(false);
      
      // Clear the file input
      if (fileInputRef.current) {
        fileInputRef.current.value = '';
      }
      
      // Call the callback if provided
      if (onImportComplete) {
        onImportComplete(importedTemplateId);
      }
    } catch (error) {
      handleError(error, 'Failed to import template');
      setImportError('Failed to import template. Please try again.');
    } finally {
      setIsImporting(false);
    }
  };
  
  // Handle template export
  const handleExport = async () => {
    if (!templateId || !template) {
      toast({
        title: 'Export Failed',
        description: 'No template selected for export.',
        variant: 'destructive'
      });
      return;
    }
    
    setIsExporting(true);
    
    try {
      // Get the template data for export
      const exportData = await exportTemplateMutation({ templateId });
      
      // Create a Blob with the JSON data
      const blob = new Blob([JSON.stringify(exportData, null, 2)], { type: 'application/json' });
      
      // Create a download link and trigger it
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = `protocol-template-${template.name.toLowerCase().replace(/\s+/g, '-')}-v${template.version}.json`;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      
      toast({
        title: 'Template Exported',
        description: 'Protocol template has been exported successfully.',
      });
      
      // Call the callback if provided
      if (onExportComplete) {
        onExportComplete();
      }
    } catch (error) {
      handleError(error, 'Failed to export template');
      toast({
        title: 'Export Failed',
        description: 'Failed to export template. Please try again.',
        variant: 'destructive'
      });
    } finally {
      setIsExporting(false);
    }
  };
  
  return (
    <div className={`space-y-6 ${className}`}>
      {/* Import Section */}
      {showImport && (
        <div className="space-y-4">
          <h3 className="text-lg font-medium">Import Protocol Template</h3>
          
          <div className="flex flex-col space-y-4">
            <div className="flex items-center space-x-4">
              <input
                type="file"
                ref={fileInputRef}
                accept=".json"
                onChange={handleFileSelect}
                className="hidden"
                id="template-file-input"
              />
              <button
                onClick={() => fileInputRef.current?.click()}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
                disabled={isImporting}
              >
                Select Template File
              </button>
              
              {importFile && (
                <span className="text-sm text-gray-600">
                  {importFile.name} ({(importFile.size / 1024).toFixed(1)} KB)
                </span>
              )}
            </div>
            
            {importError && (
              <div className="px-4 py-3 bg-red-100 text-red-800 rounded-md text-sm">
                {importError}
              </div>
            )}
          </div>
          
          {/* Import Preview */}
          {showImportPreview && importData && (
            <div className="border border-gray-300 rounded-md p-4 space-y-4 bg-gray-50">
              <h4 className="font-medium">Import Preview</h4>
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div>
                  <p className="text-sm font-medium text-gray-600">Template Name</p>
                  <p>{importData.template.name}</p>
                </div>
                
                <div>
                  <p className="text-sm font-medium text-gray-600">Version</p>
                  <p>{importData.template.version}</p>
                </div>
                
                {importData.template.description && (
                  <div className="col-span-1 md:col-span-2">
                    <p className="text-sm font-medium text-gray-600">Description</p>
                    <p>{importData.template.description}</p>
                  </div>
                )}
                
                <div>
                  <p className="text-sm font-medium text-gray-600">Steps</p>
                  <p>{importData.template.steps.length} steps</p>
                </div>
                
                <div>
                  <p className="text-sm font-medium text-gray-600">Export Date</p>
                  <p>{new Date(importData.exportDate).toLocaleString()}</p>
                </div>
                
                {importData.template.tags && importData.template.tags.length > 0 && (
                  <div className="col-span-1 md:col-span-2">
                    <p className="text-sm font-medium text-gray-600">Tags</p>
                    <div className="flex flex-wrap gap-2 mt-1">
                      {importData.template.tags.map((tag, index) => (
                        <span 
                          key={index}
                          className="px-2 py-1 bg-blue-100 text-blue-800 rounded-full text-xs"
                        >
                          {tag}
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
              
              <div className="flex justify-end space-x-4 pt-2">
                <button
                  onClick={() => {
                    setImportFile(null);
                    setImportData(null);
                    setShowImportPreview(false);
                    if (fileInputRef.current) {
                      fileInputRef.current.value = '';
                    }
                  }}
                  className="px-4 py-2 border border-gray-300 text-gray-700 rounded-md hover:bg-gray-100 transition-colors"
                  disabled={isImporting}
                >
                  Cancel
                </button>
                <button
                  onClick={handleImport}
                  className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
                  disabled={isImporting}
                >
                  {isImporting ? 'Importing...' : 'Import Template'}
                </button>
              </div>
            </div>
          )}
        </div>
      )}
      
      {/* Export Section */}
      {showExport && templateId && template && (
        <div className="space-y-4">
          <h3 className="text-lg font-medium">Export Protocol Template</h3>
          
          <div className="border border-gray-300 rounded-md p-4 space-y-4 bg-gray-50">
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div>
                <p className="text-sm font-medium text-gray-600">Template Name</p>
                <p>{template.name}</p>
              </div>
              
              <div>
                <p className="text-sm font-medium text-gray-600">Version</p>
                <p>{template.version}</p>
              </div>
              
              {template.description && (
                <div className="col-span-1 md:col-span-2">
                  <p className="text-sm font-medium text-gray-600">Description</p>
                  <p>{template.description}</p>
                </div>
              )}
              
              <div>
                <p className="text-sm font-medium text-gray-600">Steps</p>
                <p>{template.steps.length} steps</p>
              </div>
              
              {template.tags && template.tags.length > 0 && (
                <div className="col-span-1 md:col-span-2">
                  <p className="text-sm font-medium text-gray-600">Tags</p>
                  <div className="flex flex-wrap gap-2 mt-1">
                    {template.tags.map((tag, index) => (
                      <span 
                        key={index}
                        className="px-2 py-1 bg-blue-100 text-blue-800 rounded-full text-xs"
                      >
                        {tag}
                      </span>
                    ))}
                  </div>
                </div>
              )}
            </div>
            
            <div className="flex justify-end pt-2">
              <button
                onClick={handleExport}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
                disabled={isExporting}
              >
                {isExporting ? 'Exporting...' : 'Export Template'}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}