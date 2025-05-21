/**
 * Protocol Template Management Component
 * Provides UI for creating, managing, and using protocol templates
 */

import React, { useState, useEffect } from 'react';
import { useRouter } from 'next/router';
import { useProtocolTemplates, TemplateCategory, useTemplateVersions } from '../hooks/use-convex-protocol-templates';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../ui/tabs';
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from '../../ui/card';
import { Button } from '../../ui/button';
import { Badge } from '../../ui/badge';
import { Input } from '../../ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '../../ui/select';
import { PlusIcon, SearchIcon, FilterIcon, InfoIcon, BookmarkIcon, ClockIcon, PencilIcon, TrashIcon, CopyIcon } from '../../ui/icons';
import { Dialog, DialogContent, DialogDescription, DialogFooter, DialogHeader, DialogTitle, DialogTrigger } from '../../ui/dialog';
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from '../../ui/tooltip';
import { Alert, AlertDescription, AlertTitle } from '../../ui/alert';

export interface TemplateManagementProps {
  userId?: string;
  userName?: string;
}

export function TemplateManagement({ userId, userName }: TemplateManagementProps) {
  const router = useRouter();
  const [activeTab, setActiveTab] = useState('all');
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedCategory, setSelectedCategory] = useState<TemplateCategory | ''>('');
  const [showDialog, setShowDialog] = useState(false);
  const [selectedTemplate, setSelectedTemplate] = useState<any>(null);
  const [customizations, setCustomizations] = useState({
    name: '',
    description: '',
  });

  // Fetch templates with filtering options
  const templateFilter = {
    ...(searchQuery ? { name: searchQuery } : {}),
    ...(selectedCategory ? { category: selectedCategory as TemplateCategory } : {}),
  };
  
  const options = {
    limit: 50,
    sortBy: 'updatedAt',
    sortDirection: 'desc' as const,
  };
  
  const { 
    templates, 
    loading, 
    createFromTemplate, 
    getTemplateCategories 
  } = useProtocolTemplates(templateFilter, options);

  // Filter templates based on the active tab
  const filteredTemplates = templates.filter(template => {
    if (activeTab === 'my' && template.createdBy?.toString() !== userId) {
      return false;
    }
    if (activeTab === 'public') {
      return template.public === true;
    }
    return true;
  });

  // Fetch template versions if a template is selected
  const { versions } = useTemplateVersions(selectedTemplate?.name);

  // Handle creating a protocol from template
  const handleCreateFromTemplate = async () => {
    if (!selectedTemplate || !customizations.name) return;
    
    try {
      const protocolId = await createFromTemplate(selectedTemplate._id, {
        name: customizations.name,
        description: customizations.description,
      });
      
      setShowDialog(false);
      router.push(`/protocols/${protocolId}`);
    } catch (error) {
      console.error('Error creating protocol from template:', error);
    }
  };

  // Display specific template information
  const openTemplateDetails = (template: any) => {
    setSelectedTemplate(template);
    setCustomizations({
      name: `${template.name} - My Version`,
      description: template.description || '',
    });
    setShowDialog(true);
  };

  // Helper function to estimate protocol duration
  const estimateProtocolDuration = (steps: any[]) => {
    const totalMinutes = steps.reduce((total, step) => {
      if (step.duration) {
        if (step.durationUnit === 'hours') {
          return total + (step.duration * 60);
        } else if (step.durationUnit === 'days') {
          return total + (step.duration * 24 * 60);
        }
        return total + step.duration;
      }
      return total;
    }, 0);
    
    if (totalMinutes >= 24 * 60) {
      return `${Math.round(totalMinutes / (24 * 60))} days`;
    } else if (totalMinutes >= 60) {
      return `${Math.round(totalMinutes / 60)} hours`;
    }
    return `${totalMinutes} min`;
  };

  return (
    <div className="space-y-6">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center">
        <div>
          <h1 className="text-3xl font-bold mb-2">Protocol Templates</h1>
          <p className="text-muted-foreground">
            Standard templates for reproducible cryopreservation protocols
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Button 
            onClick={() => router.push('/protocols/templates/create')}
            className="inline-flex items-center justify-center gap-2"
          >
            <PlusIcon className="w-4 h-4" />
            Create Template
          </Button>
        </div>
      </div>
      
      {/* Filter and search bar */}
      <div className="flex flex-col md:flex-row gap-4 items-center">
        <div className="relative flex-1">
          <SearchIcon className="absolute left-3 top-1/2 transform -translate-y-1/2 text-muted-foreground w-4 h-4" />
          <Input
            type="search"
            placeholder="Search templates..."
            className="pl-10"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
          />
        </div>
        
        <div className="w-full md:w-auto">
          <Select
            value={selectedCategory}
            onValueChange={(value) => setSelectedCategory(value as TemplateCategory | '')}
          >
            <SelectTrigger className="w-full md:w-[180px]">
              <div className="flex items-center gap-2">
                <FilterIcon className="w-4 h-4" />
                <SelectValue placeholder="Category" />
              </div>
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="">All Categories</SelectItem>
              {getTemplateCategories().map(category => (
                <SelectItem key={category} value={category}>{category}</SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>
      </div>
      
      {/* Tabs */}
      <Tabs defaultValue="all" value={activeTab} onValueChange={setActiveTab}>
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="all">All Templates</TabsTrigger>
          <TabsTrigger value="my">My Templates</TabsTrigger>
          <TabsTrigger value="public">Public Templates</TabsTrigger>
        </TabsList>
        
        <TabsContent value={activeTab} className="mt-6">
          {loading ? (
            <div className="py-12 text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto"></div>
              <p className="mt-4 text-muted-foreground">Loading templates...</p>
            </div>
          ) : filteredTemplates.length === 0 ? (
            <div className="py-12 text-center border rounded-lg">
              <p className="text-muted-foreground mb-2">No templates found</p>
              <Button 
                variant="outline" 
                size="sm"
                onClick={() => {
                  setSearchQuery('');
                  setSelectedCategory('');
                }}
              >
                Clear filters
              </Button>
            </div>
          ) : (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {filteredTemplates.map((template) => (
                <Card key={template._id} className="overflow-hidden">
                  <CardHeader className="pb-3">
                    <div className="flex items-center justify-between">
                      <CardTitle className="text-lg font-semibold">{template.name}</CardTitle>
                      {template.public && (
                        <TooltipProvider>
                          <Tooltip>
                            <TooltipTrigger>
                              <Badge variant="outline" className="ml-2">
                                Public
                              </Badge>
                            </TooltipTrigger>
                            <TooltipContent>
                              <p>This template is available to all users</p>
                            </TooltipContent>
                          </Tooltip>
                        </TooltipProvider>
                      )}
                    </div>
                    <div className="flex flex-wrap items-center gap-2 text-xs text-muted-foreground">
                      <Badge variant="secondary" className="rounded-full">
                        v{template.version}
                      </Badge>
                      {template.category && (
                        <Badge variant="outline" className="rounded-full">
                          {template.category}
                        </Badge>
                      )}
                    </div>
                    <CardDescription className="mt-2 line-clamp-2">
                      {template.description || "No description provided"}
                    </CardDescription>
                  </CardHeader>
                  
                  <CardContent className="text-sm text-muted-foreground pb-3">
                    <div className="flex flex-wrap gap-4">
                      <div className="flex items-center gap-1">
                        <ClockIcon className="w-4 h-4" />
                        {estimateProtocolDuration(template.steps)}
                      </div>
                      <div className="flex items-center gap-1">
                        <InfoIcon className="w-4 h-4" />
                        {template.steps.length} steps
                      </div>
                    </div>
                  </CardContent>
                  
                  <CardFooter className="flex justify-between pt-3 border-t">
                    <div className="text-xs text-muted-foreground">
                      Updated {new Date(template.updatedAt).toLocaleDateString()}
                    </div>
                    <div className="flex gap-2">
                      <Button 
                        variant="ghost" 
                        size="sm"
                        onClick={() => router.push(`/protocols/templates/${template._id}`)}
                      >
                        <PencilIcon className="w-4 h-4 mr-1" />
                        Edit
                      </Button>
                      <Button 
                        variant="default" 
                        size="sm"
                        onClick={() => openTemplateDetails(template)}
                      >
                        <CopyIcon className="w-4 h-4 mr-1" />
                        Use
                      </Button>
                    </div>
                  </CardFooter>
                </Card>
              ))}
            </div>
          )}
        </TabsContent>
      </Tabs>
      
      {/* Dialog for template details and customization */}
      {selectedTemplate && (
        <Dialog open={showDialog} onOpenChange={setShowDialog}>
          <DialogContent className="sm:max-w-lg">
            <DialogHeader>
              <DialogTitle>Create Protocol from Template</DialogTitle>
              <DialogDescription>
                Customize this template to create your own protocol
              </DialogDescription>
            </DialogHeader>
            
            <div className="grid gap-4 py-4">
              <div>
                <h3 className="font-medium mb-2">Template Details</h3>
                <Card className="bg-muted/40">
                  <CardContent className="pt-6">
                    <div className="space-y-2">
                      <div className="flex justify-between">
                        <div className="font-medium">{selectedTemplate.name}</div>
                        <Badge variant="outline">v{selectedTemplate.version}</Badge>
                      </div>
                      <p className="text-sm text-muted-foreground">
                        {selectedTemplate.description}
                      </p>
                      <div className="text-sm flex gap-4">
                        <div>{selectedTemplate.steps.length} steps</div>
                        <div>{estimateProtocolDuration(selectedTemplate.steps)}</div>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              </div>
              
              {versions && versions.length > 1 && (
                <Alert>
                  <BookmarkIcon className="h-4 w-4" />
                  <AlertTitle>Multiple versions available</AlertTitle>
                  <AlertDescription>
                    This template has {versions.length} versions. You are using v{selectedTemplate.version}.
                  </AlertDescription>
                </Alert>
              )}
              
              <div className="space-y-2">
                <label htmlFor="protocol-name" className="text-sm font-medium">
                  Protocol Name <span className="text-red-500">*</span>
                </label>
                <Input
                  id="protocol-name"
                  value={customizations.name}
                  onChange={(e) => setCustomizations({...customizations, name: e.target.value})}
                  placeholder="Enter protocol name"
                />
              </div>
              
              <div className="space-y-2">
                <label htmlFor="protocol-description" className="text-sm font-medium">
                  Description
                </label>
                <textarea
                  id="protocol-description"
                  className="min-h-24 w-full rounded-md border border-input bg-background px-3 py-2 text-sm"
                  value={customizations.description}
                  onChange={(e) => setCustomizations({...customizations, description: e.target.value})}
                  placeholder="Enter protocol description"
                />
              </div>
            </div>
            
            <DialogFooter>
              <Button variant="outline" onClick={() => setShowDialog(false)}>
                Cancel
              </Button>
              <Button 
                onClick={handleCreateFromTemplate}
                disabled={!customizations.name}
              >
                Create Protocol
              </Button>
            </DialogFooter>
          </DialogContent>
        </Dialog>
      )}
    </div>
  );
}

export default TemplateManagement;