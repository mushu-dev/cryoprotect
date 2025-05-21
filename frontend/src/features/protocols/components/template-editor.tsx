/**
 * Protocol Template Editor Component
 * Provides UI for creating and editing protocol templates
 */

import React, { useState, useEffect } from 'react';
import { useRouter } from 'next/router';
import { useProtocolTemplates, TemplateCategory } from '../hooks/use-convex-protocol-templates';
import { Button } from '../../ui/button';
import { Input } from '../../ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '../../ui/select';
import { Card, CardContent, CardHeader, CardTitle, CardDescription, CardFooter } from '../../ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../ui/tabs';
import { Switch } from '../../ui/switch';
import { Label } from '../../ui/label';
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from '../../ui/accordion';
import { ChevronLeftIcon, PlusIcon, TrashIcon, CopyIcon, ArrowUpIcon, ArrowDownIcon } from '../../ui/icons';
import { useToast } from '../../ui/use-toast';
import { Dialog, DialogContent, DialogDescription, DialogFooter, DialogHeader, DialogTitle, DialogTrigger } from '../../ui/dialog';
import { Alert, AlertDescription, AlertTitle } from '../../ui/alert';
import { Id } from "../../../../convex/_generated/dataModel";

export interface TemplateEditorProps {
  templateId?: string;
  userId?: string;
  userName?: string;
}

// Step interface matching backend schema
interface ProtocolStep {
  id?: string;
  name: string;
  description?: string;
  parameters?: Record<string, any>;
  duration?: number;
  durationUnit?: string;
  temperature?: number;
  temperatureUnit?: string;
}

export function TemplateEditor({ templateId, userId, userName }: TemplateEditorProps) {
  const router = useRouter();
  const { toast } = useToast();
  const { getTemplate, createTemplate, updateTemplate, getTemplateCategories } = useProtocolTemplates();
  
  const [loading, setLoading] = useState(templateId ? true : false);
  const [showConfirmDialog, setShowConfirmDialog] = useState(false);
  const [activeTab, setActiveTab] = useState('basic');
  
  // Template state
  const [template, setTemplate] = useState<{
    name: string;
    description: string;
    steps: ProtocolStep[];
    parameters: Record<string, any>;
    category: TemplateCategory | '';
    public: boolean;
    version: number;
  }>({
    name: '',
    description: '',
    steps: [],
    parameters: {},
    category: '',
    public: false,
    version: 1
  });
  
  // Step editor state
  const [editingStepIndex, setEditingStepIndex] = useState<number | null>(null);
  const [stepEditor, setStepEditor] = useState<ProtocolStep>({
    name: '',
    description: '',
    duration: 0,
    durationUnit: 'minutes',
    temperature: 20,
    temperatureUnit: 'C',
    parameters: {}
  });
  
  // Parameter editor state
  const [parameterName, setParameterName] = useState('');
  const [parameterValue, setParameterValue] = useState('');
  
  // Load template if editing
  useEffect(() => {
    if (templateId) {
      const loadTemplate = async () => {
        try {
          setLoading(true);
          const templateData = await getTemplate({ protocolId: templateId as Id<"protocols"> });
          
          if (templateData) {
            setTemplate({
              name: templateData.name,
              description: templateData.description || '',
              steps: templateData.steps || [],
              parameters: templateData.parameters || {},
              category: templateData.category || '',
              public: templateData.public || false,
              version: templateData.version || 1
            });
          }
        } catch (error) {
          console.error('Error loading template:', error);
          toast({
            title: 'Error',
            description: 'Failed to load template',
            variant: 'destructive'
          });
        } finally {
          setLoading(false);
        }
      };
      
      loadTemplate();
    }
  }, [templateId, getTemplate, toast]);
  
  // Handle save/update
  const handleSave = async () => {
    try {
      if (!template.name) {
        toast({
          title: 'Validation Error',
          description: 'Template name is required',
          variant: 'destructive'
        });
        return;
      }
      
      if (template.steps.length === 0) {
        toast({
          title: 'Validation Error',
          description: 'Template must have at least one step',
          variant: 'destructive'
        });
        return;
      }
      
      if (templateId) {
        // Update existing template
        await updateTemplate(templateId as Id<"protocols">, {
          name: template.name,
          description: template.description,
          steps: template.steps,
          parameters: template.parameters,
          category: template.category || undefined,
          public: template.public
        });
        
        toast({
          title: 'Success',
          description: 'Template updated successfully'
        });
      } else {
        // Create new template
        const newTemplateId = await createTemplate({
          name: template.name,
          description: template.description,
          steps: template.steps,
          parameters: template.parameters,
          category: template.category || undefined,
          public: template.public,
          isTemplate: true
        });
        
        toast({
          title: 'Success',
          description: 'Template created successfully'
        });
        
        // Redirect to template list
        router.push('/protocols/templates');
      }
    } catch (error) {
      console.error('Error saving template:', error);
      toast({
        title: 'Error',
        description: 'Failed to save template',
        variant: 'destructive'
      });
    }
  };
  
  // Step management
  const openStepEditor = (index: number) => {
    setEditingStepIndex(index);
    setStepEditor(index < 0 ? {
      name: '',
      description: '',
      duration: 0,
      durationUnit: 'minutes',
      temperature: 20,
      temperatureUnit: 'C',
      parameters: {}
    } : { ...template.steps[index] });
  };
  
  const saveStep = () => {
    if (!stepEditor.name) {
      toast({
        title: 'Validation Error',
        description: 'Step name is required',
        variant: 'destructive'
      });
      return;
    }
    
    const updatedSteps = [...template.steps];
    
    if (editingStepIndex !== null && editingStepIndex >= 0) {
      // Update existing step
      updatedSteps[editingStepIndex] = { ...stepEditor };
    } else {
      // Add new step
      updatedSteps.push({ ...stepEditor });
    }
    
    setTemplate({
      ...template,
      steps: updatedSteps
    });
    
    setEditingStepIndex(null);
  };
  
  const removeStep = (index: number) => {
    const updatedSteps = [...template.steps];
    updatedSteps.splice(index, 1);
    
    setTemplate({
      ...template,
      steps: updatedSteps
    });
  };
  
  const moveStep = (index: number, direction: 'up' | 'down') => {
    if (
      (direction === 'up' && index === 0) || 
      (direction === 'down' && index === template.steps.length - 1)
    ) {
      return;
    }
    
    const updatedSteps = [...template.steps];
    const targetIndex = direction === 'up' ? index - 1 : index + 1;
    
    [updatedSteps[index], updatedSteps[targetIndex]] = [updatedSteps[targetIndex], updatedSteps[index]];
    
    setTemplate({
      ...template,
      steps: updatedSteps
    });
  };
  
  // Parameter management
  const addParameter = () => {
    if (!parameterName) {
      toast({
        title: 'Validation Error',
        description: 'Parameter name is required',
        variant: 'destructive'
      });
      return;
    }
    
    setTemplate({
      ...template,
      parameters: {
        ...template.parameters,
        [parameterName]: parameterValue
      }
    });
    
    setParameterName('');
    setParameterValue('');
  };
  
  const removeParameter = (key: string) => {
    const updatedParameters = { ...template.parameters };
    delete updatedParameters[key];
    
    setTemplate({
      ...template,
      parameters: updatedParameters
    });
  };
  
  // Step parameter management
  const addStepParameter = (name: string, value: string) => {
    if (!name) return;
    
    setStepEditor({
      ...stepEditor,
      parameters: {
        ...(stepEditor.parameters || {}),
        [name]: value
      }
    });
  };
  
  const removeStepParameter = (key: string) => {
    const updatedParameters = { ...(stepEditor.parameters || {}) };
    delete updatedParameters[key];
    
    setStepEditor({
      ...stepEditor,
      parameters: updatedParameters
    });
  };
  
  // Handle cancel
  const handleCancel = () => {
    if (
      template.name || 
      template.description || 
      template.steps.length > 0 || 
      Object.keys(template.parameters).length > 0
    ) {
      setShowConfirmDialog(true);
    } else {
      router.push('/protocols/templates');
    }
  };
  
  if (loading) {
    return (
      <div className="py-8 text-center">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto"></div>
        <p className="mt-4 text-muted-foreground">Loading template...</p>
      </div>
    );
  }
  
  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div className="flex items-center space-x-2">
          <Button 
            variant="ghost" 
            size="icon"
            onClick={handleCancel}
          >
            <ChevronLeftIcon className="h-5 w-5" />
            <span className="sr-only">Back</span>
          </Button>
          <h1 className="text-3xl font-bold">
            {templateId ? 'Edit Template' : 'Create Template'}
          </h1>
        </div>
        
        <div className="flex space-x-2">
          <Button variant="outline" onClick={handleCancel}>
            Cancel
          </Button>
          <Button onClick={handleSave}>
            {templateId ? 'Update' : 'Save'} Template
          </Button>
        </div>
      </div>
      
      <Tabs value={activeTab} onValueChange={setActiveTab}>
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="basic">Basic Information</TabsTrigger>
          <TabsTrigger value="steps">Protocol Steps</TabsTrigger>
          <TabsTrigger value="parameters">Parameters</TabsTrigger>
        </TabsList>
        
        {/* Basic Information Tab */}
        <TabsContent value="basic" className="space-y-6 py-4">
          <div className="grid gap-6">
            <div className="space-y-2">
              <Label htmlFor="template-name">
                Template Name <span className="text-red-500">*</span>
              </Label>
              <Input
                id="template-name"
                value={template.name}
                onChange={(e) => setTemplate({...template, name: e.target.value})}
                placeholder="Enter template name"
              />
            </div>
            
            <div className="space-y-2">
              <Label htmlFor="template-description">Description</Label>
              <textarea
                id="template-description"
                className="min-h-32 w-full rounded-md border border-input bg-background px-3 py-2 text-sm"
                value={template.description}
                onChange={(e) => setTemplate({...template, description: e.target.value})}
                placeholder="Enter template description"
              />
            </div>
            
            <div className="space-y-2">
              <Label htmlFor="template-category">Category</Label>
              <Select
                value={template.category}
                onValueChange={(value) => setTemplate({...template, category: value as TemplateCategory | ''})}
              >
                <SelectTrigger id="template-category">
                  <SelectValue placeholder="Select category" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="">None</SelectItem>
                  {getTemplateCategories().map(category => (
                    <SelectItem key={category} value={category}>{category}</SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>
            
            <div className="flex items-center space-x-2">
              <Switch
                id="template-public"
                checked={template.public}
                onCheckedChange={(checked) => setTemplate({...template, public: checked})}
              />
              <Label htmlFor="template-public">Make template public</Label>
            </div>
          </div>
        </TabsContent>
        
        {/* Protocol Steps Tab */}
        <TabsContent value="steps" className="space-y-6 py-4">
          {template.steps.length === 0 ? (
            <div className="py-12 text-center border rounded-lg">
              <p className="text-muted-foreground mb-2">No steps added yet</p>
              <Button onClick={() => openStepEditor(-1)}>
                <PlusIcon className="h-4 w-4 mr-2" />
                Add First Step
              </Button>
            </div>
          ) : (
            <div className="space-y-4">
              <Accordion type="multiple" className="w-full">
                {template.steps.map((step, index) => (
                  <AccordionItem key={index} value={`step-${index}`}>
                    <AccordionTrigger className="text-left">
                      <div className="flex items-center">
                        <span className="bg-muted rounded-full w-6 h-6 inline-flex items-center justify-center mr-2">
                          {index + 1}
                        </span>
                        {step.name}
                      </div>
                    </AccordionTrigger>
                    <AccordionContent>
                      <div className="p-4 space-y-4">
                        {step.description && (
                          <p className="text-sm text-muted-foreground">{step.description}</p>
                        )}
                        
                        <div className="grid grid-cols-2 gap-4 text-sm">
                          {step.duration && (
                            <div>
                              <span className="font-medium">Duration:</span> {step.duration} {step.durationUnit}
                            </div>
                          )}
                          
                          {step.temperature && (
                            <div>
                              <span className="font-medium">Temperature:</span> {step.temperature}°{step.temperatureUnit}
                            </div>
                          )}
                        </div>
                        
                        {step.parameters && Object.keys(step.parameters).length > 0 && (
                          <div className="space-y-2">
                            <h4 className="text-sm font-medium">Parameters:</h4>
                            <div className="bg-muted rounded-md p-2">
                              {Object.entries(step.parameters).map(([key, value]) => (
                                <div key={key} className="text-sm">
                                  <span className="font-medium">{key}:</span> {JSON.stringify(value)}
                                </div>
                              ))}
                            </div>
                          </div>
                        )}
                        
                        <div className="flex justify-end space-x-2 pt-2">
                          <Button 
                            variant="outline" 
                            size="sm"
                            onClick={() => moveStep(index, 'up')}
                            disabled={index === 0}
                          >
                            <ArrowUpIcon className="h-4 w-4" />
                            <span className="sr-only">Move Up</span>
                          </Button>
                          <Button 
                            variant="outline" 
                            size="sm"
                            onClick={() => moveStep(index, 'down')}
                            disabled={index === template.steps.length - 1}
                          >
                            <ArrowDownIcon className="h-4 w-4" />
                            <span className="sr-only">Move Down</span>
                          </Button>
                          <Button 
                            variant="outline" 
                            size="sm"
                            onClick={() => openStepEditor(index)}
                          >
                            Edit
                          </Button>
                          <Button
                            variant="destructive"
                            size="sm"
                            onClick={() => removeStep(index)}
                          >
                            <TrashIcon className="h-4 w-4" />
                            <span className="sr-only">Delete</span>
                          </Button>
                        </div>
                      </div>
                    </AccordionContent>
                  </AccordionItem>
                ))}
              </Accordion>
              
              <div className="flex justify-center pt-4">
                <Button onClick={() => openStepEditor(-1)}>
                  <PlusIcon className="h-4 w-4 mr-2" />
                  Add Step
                </Button>
              </div>
            </div>
          )}
        </TabsContent>
        
        {/* Parameters Tab */}
        <TabsContent value="parameters" className="space-y-6 py-4">
          <div className="space-y-6">
            <div className="flex gap-2 items-end">
              <div className="flex-1 space-y-2">
                <Label htmlFor="parameter-name">Parameter Name</Label>
                <Input
                  id="parameter-name"
                  value={parameterName}
                  onChange={(e) => setParameterName(e.target.value)}
                  placeholder="Enter parameter name"
                />
              </div>
              <div className="flex-1 space-y-2">
                <Label htmlFor="parameter-value">Value</Label>
                <Input
                  id="parameter-value"
                  value={parameterValue}
                  onChange={(e) => setParameterValue(e.target.value)}
                  placeholder="Enter parameter value"
                />
              </div>
              <Button onClick={addParameter}>
                Add
              </Button>
            </div>
            
            {Object.keys(template.parameters).length === 0 ? (
              <div className="py-8 text-center border rounded-lg">
                <p className="text-muted-foreground">No parameters added yet</p>
              </div>
            ) : (
              <Card>
                <CardHeader>
                  <CardTitle>Template Parameters</CardTitle>
                  <CardDescription>
                    Global parameters for the entire protocol
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="space-y-2">
                    {Object.entries(template.parameters).map(([key, value]) => (
                      <div key={key} className="flex justify-between items-center p-2 border-b last:border-0">
                        <div>
                          <span className="font-medium">{key}:</span> {JSON.stringify(value)}
                        </div>
                        <Button 
                          variant="ghost" 
                          size="sm"
                          onClick={() => removeParameter(key)}
                        >
                          <TrashIcon className="h-4 w-4 text-red-500" />
                          <span className="sr-only">Remove</span>
                        </Button>
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
        </TabsContent>
      </Tabs>
      
      {/* Step Editor Dialog */}
      {editingStepIndex !== null && (
        <Dialog open={editingStepIndex !== null} onOpenChange={() => setEditingStepIndex(null)}>
          <DialogContent className="sm:max-w-lg">
            <DialogHeader>
              <DialogTitle>
                {editingStepIndex >= 0 ? 'Edit Step' : 'Add Step'}
              </DialogTitle>
              <DialogDescription>
                Define the parameters for this protocol step
              </DialogDescription>
            </DialogHeader>
            
            <div className="grid gap-4 py-4">
              <div className="space-y-2">
                <Label htmlFor="step-name">
                  Step Name <span className="text-red-500">*</span>
                </Label>
                <Input
                  id="step-name"
                  value={stepEditor.name}
                  onChange={(e) => setStepEditor({...stepEditor, name: e.target.value})}
                  placeholder="Enter step name"
                />
              </div>
              
              <div className="space-y-2">
                <Label htmlFor="step-description">Description</Label>
                <textarea
                  id="step-description"
                  className="min-h-24 w-full rounded-md border border-input bg-background px-3 py-2 text-sm"
                  value={stepEditor.description || ''}
                  onChange={(e) => setStepEditor({...stepEditor, description: e.target.value})}
                  placeholder="Enter step description"
                />
              </div>
              
              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="step-duration">Duration</Label>
                  <div className="flex space-x-2">
                    <Input
                      id="step-duration"
                      type="number"
                      value={stepEditor.duration || ''}
                      onChange={(e) => setStepEditor({
                        ...stepEditor, 
                        duration: e.target.value ? Number(e.target.value) : undefined
                      })}
                      placeholder="Duration"
                    />
                    <Select
                      value={stepEditor.durationUnit || 'minutes'}
                      onValueChange={(value) => setStepEditor({...stepEditor, durationUnit: value})}
                    >
                      <SelectTrigger className="w-[110px]">
                        <SelectValue placeholder="Unit" />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="seconds">seconds</SelectItem>
                        <SelectItem value="minutes">minutes</SelectItem>
                        <SelectItem value="hours">hours</SelectItem>
                        <SelectItem value="days">days</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>
                
                <div className="space-y-2">
                  <Label htmlFor="step-temperature">Temperature</Label>
                  <div className="flex space-x-2">
                    <Input
                      id="step-temperature"
                      type="number"
                      value={stepEditor.temperature || ''}
                      onChange={(e) => setStepEditor({
                        ...stepEditor, 
                        temperature: e.target.value ? Number(e.target.value) : undefined
                      })}
                      placeholder="Temperature"
                    />
                    <Select
                      value={stepEditor.temperatureUnit || 'C'}
                      onValueChange={(value) => setStepEditor({...stepEditor, temperatureUnit: value})}
                    >
                      <SelectTrigger className="w-[60px]">
                        <SelectValue placeholder="Unit" />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="C">°C</SelectItem>
                        <SelectItem value="F">°F</SelectItem>
                        <SelectItem value="K">K</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>
              </div>
              
              <div className="space-y-4">
                <div className="flex justify-between items-center">
                  <Label>Step Parameters</Label>
                </div>
                
                <div className="flex gap-2 items-end">
                  <div className="flex-1">
                    <Input
                      placeholder="Parameter name"
                      value={parameterName}
                      onChange={(e) => setParameterName(e.target.value)}
                    />
                  </div>
                  <div className="flex-1">
                    <Input
                      placeholder="Value"
                      value={parameterValue}
                      onChange={(e) => setParameterValue(e.target.value)}
                    />
                  </div>
                  <Button 
                    size="sm"
                    onClick={() => {
                      addStepParameter(parameterName, parameterValue);
                      setParameterName('');
                      setParameterValue('');
                    }}
                  >
                    Add
                  </Button>
                </div>
                
                {stepEditor.parameters && Object.keys(stepEditor.parameters).length > 0 ? (
                  <div className="border rounded-md divide-y">
                    {Object.entries(stepEditor.parameters).map(([key, value]) => (
                      <div key={key} className="flex justify-between items-center p-2">
                        <div>
                          <span className="font-medium">{key}:</span> {JSON.stringify(value)}
                        </div>
                        <Button 
                          variant="ghost" 
                          size="sm" 
                          onClick={() => removeStepParameter(key)}
                        >
                          <TrashIcon className="h-4 w-4 text-red-500" />
                        </Button>
                      </div>
                    ))}
                  </div>
                ) : (
                  <div className="text-center py-2 text-sm text-muted-foreground border rounded-md">
                    No parameters added
                  </div>
                )}
              </div>
            </div>
            
            <DialogFooter>
              <Button variant="outline" onClick={() => setEditingStepIndex(null)}>
                Cancel
              </Button>
              <Button onClick={saveStep}>
                {editingStepIndex >= 0 ? 'Update' : 'Add'} Step
              </Button>
            </DialogFooter>
          </DialogContent>
        </Dialog>
      )}
      
      {/* Confirmation Dialog */}
      <Dialog open={showConfirmDialog} onOpenChange={setShowConfirmDialog}>
        <DialogContent>
          <DialogHeader>
            <DialogTitle>Discard changes?</DialogTitle>
            <DialogDescription>
              Any unsaved changes will be lost. Are you sure you want to leave?
            </DialogDescription>
          </DialogHeader>
          <DialogFooter>
            <Button variant="outline" onClick={() => setShowConfirmDialog(false)}>
              Cancel
            </Button>
            <Button 
              variant="destructive"
              onClick={() => {
                setShowConfirmDialog(false);
                router.push('/protocols/templates');
              }}
            >
              Discard Changes
            </Button>
          </DialogFooter>
        </DialogContent>
      </Dialog>
    </div>
  );
}

export default TemplateEditor;