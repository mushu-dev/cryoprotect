'use client'

import { useState } from 'react'
import { useForm, useFieldArray } from 'react-hook-form'
import { useRouter } from 'next/navigation'
import { useCreateMixture } from '@/features/mixtures/hooks/use-mixtures'
import { type CreateMixtureData } from '@/features/mixtures/services/mixture-service'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { Textarea } from '@/components/ui/textarea'
import { MixtureCompositionChart } from '@/features/mixtures/components/mixture-composition-chart'
import { ArrowLeft, Plus, Trash, Check } from 'lucide-react'
import Link from 'next/link'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'

interface MixtureFormComponent {
  molecule_id: string
  molecule_name?: string
  concentration: number
  concentration_unit: string
}

interface MixtureFormValues {
  name: string
  description?: string
  components: MixtureFormComponent[]
}

export default function CreateMixturePage() {
  const router = useRouter()
  const [error, setError] = useState<string | null>(null)
  
  // Mutation for creating a mixture
  const { mutate: createMixture, isPending } = useCreateMixture()
  
  // Form setup
  const form = useForm<MixtureFormValues>({
    defaultValues: {
      name: '',
      description: '',
      components: [{ 
        molecule_id: '', 
        molecule_name: 'Select a molecule',
        concentration: 100, 
        concentration_unit: '%' 
      }],
    },
  })
  
  // Setup field array for components
  const { fields, append, remove } = useFieldArray({
    control: form.control,
    name: 'components',
  })
  
  // Watch components for the chart
  const components = form.watch('components')
  
  // Chart data
  const chartComponents = components.map(component => ({
    name: component.molecule_name || 'Unnamed molecule',
    concentration: component.concentration,
    concentration_unit: component.concentration_unit,
  }))
  
  // Add a new component
  const addComponent = () => {
    // Calculate remaining concentration
    const usedConcentration = components.reduce((sum, component) => sum + (component.concentration || 0), 0)
    const remainingConcentration = Math.max(0, 100 - usedConcentration)
    
    append({ 
      molecule_id: '', 
      molecule_name: 'Select a molecule',
      concentration: remainingConcentration, 
      concentration_unit: '%' 
    })
  }
  
  // Validate the form
  const validateForm = (values: MixtureFormValues) => {
    const errors: Record<string, any> = {}
    
    // Validate name
    if (!values.name || values.name.trim() === '') {
      errors.name = 'Name is required'
    }
    
    // Validate components
    if (!values.components || values.components.length === 0) {
      errors.components = { root: 'At least one component is required' }
    } else {
      const componentErrors: Record<number, any>[] = []
      let hasErrors = false
      
      // Calculate total concentration
      const totalConcentration = values.components.reduce(
        (sum, component) => sum + (component.concentration || 0), 
        0
      )
      
      if (totalConcentration > 100) {
        errors.components = { root: 'Total concentration cannot exceed 100%' }
        hasErrors = true
      }
      
      // Check each component
      values.components.forEach((component, index) => {
        const error: Record<string, string> = {}
        
        if (!component.molecule_id) {
          error.molecule_id = 'Molecule ID is required'
          hasErrors = true
        }
        
        if (!component.concentration || component.concentration <= 0) {
          error.concentration = 'Concentration must be greater than 0'
          hasErrors = true
        }
        
        if (component.concentration > 100) {
          error.concentration = 'Concentration cannot exceed 100%'
          hasErrors = true
        }
        
        if (Object.keys(error).length > 0) {
          componentErrors[index] = error
        }
      })
      
      if (hasErrors && !errors.components) {
        errors.components = {}
      }
      
      if (componentErrors.length > 0) {
        errors.components = {
          ...errors.components,
          ...componentErrors
        }
      }
    }
    
    return errors
  }
  
  // Submit handler
  const onSubmit = (values: MixtureFormValues) => {
    setError(null)
    
    // Validate form
    const errors = validateForm(values)
    
    if (Object.keys(errors).length > 0) {
      for (const key in errors) {
        if (key === 'components' && errors.components.root) {
          setError(errors.components.root)
          return
        } else if (key === 'name') {
          setError(errors.name)
          return
        }
      }
      setError('Please fix form errors before submitting')
      return
    }
    
    // Prepare data for API
    const data: CreateMixtureData = {
      name: values.name,
      description: values.description,
      components: values.components.map(component => ({
        molecule_id: component.molecule_id,
        concentration: component.concentration,
        concentration_unit: component.concentration_unit || '%',
      })),
    }
    
    createMixture(data, {
      onSuccess: (mixture) => {
        router.push(`/mixtures/${mixture.id}`)
      },
      onError: (error) => {
        setError(error instanceof Error ? error.message : 'Failed to create mixture.')
      },
    })
  }
  
  // In a real implementation, we would use a molecule search component here.
  // For now, we'll use a simple mock with a few predefined molecules
  const mockMolecules = [
    { id: 'glycerol-id', name: 'Glycerol', formula: 'C3H8O3' },
    { id: 'dmso-id', name: 'DMSO', formula: 'C2H6OS' },
    { id: 'ethylene-glycol-id', name: 'Ethylene Glycol', formula: 'C2H6O2' },
    { id: 'propylene-glycol-id', name: 'Propylene Glycol', formula: 'C3H8O2' },
    { id: 'trehalose-id', name: 'Trehalose', formula: 'C12H22O11' },
  ]
  
  // Mock function to select a molecule - in a real implementation, this would open a search dialog
  const selectMolecule = (index: number) => {
    // Use a simple rotation through the mock molecules based on index
    const molecule = mockMolecules[index % mockMolecules.length]
    form.setValue(`components.${index}.molecule_id`, molecule.id)
    form.setValue(`components.${index}.molecule_name`, molecule.name)
  }
  
  return (
    <div className="max-w-2xl mx-auto">
      <div className="flex items-center mb-8">
        <Button 
          variant="outline" 
          size="icon" 
          className="mr-4"
          asChild
        >
          <Link href="/mixtures">
            <ArrowLeft className="h-4 w-4" />
          </Link>
        </Button>
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Create Mixture</h1>
          <p className="text-muted-foreground mt-1">
            Create a new mixture by combining multiple molecules
          </p>
        </div>
      </div>
      
      <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-8">
        {error && (
          <Alert variant="destructive">
            <AlertTitle>Error</AlertTitle>
            <AlertDescription>{error}</AlertDescription>
          </Alert>
        )}
        
        <Card>
          <CardHeader>
            <CardTitle>Mixture Details</CardTitle>
            <CardDescription>
              Enter basic information about your mixture
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="space-y-2">
              <Label htmlFor="name">Name</Label>
              <Input
                id="name"
                placeholder="Enter mixture name"
                {...form.register('name')}
              />
              {form.formState.errors.name && (
                <p className="text-destructive text-sm">
                  {form.formState.errors.name.message}
                </p>
              )}
            </div>
            
            <div className="space-y-2">
              <Label htmlFor="description">Description (optional)</Label>
              <Textarea
                id="description"
                placeholder="Enter a description of the mixture and its purpose"
                {...form.register('description')}
                rows={3}
              />
            </div>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader>
            <CardTitle>Components</CardTitle>
            <CardDescription>
              Add molecules and their concentrations to your mixture
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-6">
            {error && error.includes('concentration') && (
              <p className="text-destructive text-sm">
                {error}
              </p>
            )}
            
            {fields.map((field, index) => (
              <div key={field.id} className="p-4 border rounded-md">
                <div className="flex items-center justify-between mb-4">
                  <h3 className="font-medium">Component {index + 1}</h3>
                  {fields.length > 1 && (
                    <Button
                      type="button"
                      variant="ghost"
                      size="icon"
                      onClick={() => remove(index)}
                    >
                      <Trash className="h-4 w-4" />
                    </Button>
                  )}
                </div>
                
                <div className="space-y-4">
                  <div className="space-y-2">
                    <Label>Molecule</Label>
                    <div className="flex gap-2">
                      <Input
                        readOnly
                        value={form.watch(`components.${index}.molecule_name`)}
                        placeholder="Select a molecule"
                        className="flex-grow"
                      />
                      <Button
                        type="button"
                        variant="outline"
                        onClick={() => selectMolecule(index)}
                      >
                        Select
                      </Button>
                    </div>
                    {form.formState.errors.components?.[index]?.molecule_id && (
                      <p className="text-destructive text-sm">
                        {String(form.formState.errors.components[index]?.molecule_id?.message)}
                      </p>
                    )}
                  </div>
                  
                  <div className="space-y-2">
                    <Label htmlFor={`components.${index}.concentration`}>
                      Concentration (%)
                    </Label>
                    <Input
                      id={`components.${index}.concentration`}
                      type="number"
                      step="0.1"
                      {...form.register(`components.${index}.concentration`, { valueAsNumber: true })}
                    />
                    {form.formState.errors.components?.[index]?.concentration && (
                      <p className="text-destructive text-sm">
                        {String(form.formState.errors.components[index]?.concentration?.message)}
                      </p>
                    )}
                  </div>
                </div>
              </div>
            ))}
            
            <Button
              type="button"
              variant="outline"
              onClick={addComponent}
              className="w-full"
            >
              <Plus className="h-4 w-4 mr-2" />
              Add Component
            </Button>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader>
            <CardTitle>Mixture Preview</CardTitle>
          </CardHeader>
          <CardContent>
            {chartComponents.length > 0 ? (
              <MixtureCompositionChart 
                components={chartComponents}
                height={300}
              />
            ) : (
              <div className="flex items-center justify-center h-[300px] bg-muted/30 rounded-md">
                <p className="text-muted-foreground">Add components to preview the mixture</p>
              </div>
            )}
          </CardContent>
        </Card>
        
        <div className="flex justify-between">
          <Button variant="outline" asChild>
            <Link href="/mixtures">Cancel</Link>
          </Button>
          <Button type="submit" disabled={isPending}>
            {isPending ? 'Creating...' : (
              <>
                <Check className="mr-2 h-4 w-4" />
                Create Mixture
              </>
            )}
          </Button>
        </div>
      </form>
    </div>
  )
}