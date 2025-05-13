'use client'

import { useState } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { 
  Select, 
  SelectContent, 
  SelectItem, 
  SelectTrigger, 
  SelectValue 
} from '@/components/ui/select'
import { Search, Filter, SlidersHorizontal } from 'lucide-react'
import { Tabs, TabsList, TabsTrigger, TabsContent } from '@/components/ui/tabs'

// Mock property types
const propertyTypes = [
  'Molecular Weight',
  'LogP',
  'Hydrogen Bond Donors',
  'Hydrogen Bond Acceptors',
  'Rotatable Bonds',
  'Polar Surface Area',
  'Melting Point',
  'Boiling Point',
  'Flash Point',
  'Density',
  'Solubility',
  'Viscosity',
]

export function PropertyExplorer() {
  const [selectedProperty, setSelectedProperty] = useState<string | null>(null)
  const [searchQuery, setSearchQuery] = useState('')
  
  // Mock data for visualization - this would come from an API in a real implementation
  const mockData = {
    labels: ['Glycerol', 'DMSO', 'Ethylene Glycol', 'Propylene Glycol', 'Sucrose', 'Trehalose', 'Glucose', 'Mannose'],
    values: [92.09, 78.13, 62.07, 76.09, 342.30, 378.33, 180.16, 180.16],
    unit: 'g/mol',
  }
  
  const handleSearch = () => {
    // In a real implementation, this would trigger an API request
    console.log('Searching for:', searchQuery)
  }
  
  const handlePropertySelect = (value: string) => {
    setSelectedProperty(value)
    // In a real implementation, this would trigger an API request
  }
  
  return (
    <div className="space-y-6">
      {/* Search and filter row */}
      <div className="flex flex-col sm:flex-row gap-4">
        <div className="relative flex-grow">
          <Search className="absolute left-2.5 top-2.5 h-4 w-4 text-muted-foreground" />
          <Input 
            placeholder="Search molecules by name, formula, or property..."
            className="pl-9"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter') handleSearch()
            }}
          />
        </div>
        
        <div className="flex gap-2">
          <Button 
            variant="outline" 
            className="w-[90px]"
            onClick={handleSearch}
          >
            Search
          </Button>
          
          <Select onValueChange={handlePropertySelect}>
            <SelectTrigger className="w-[180px]">
              <Filter className="h-4 w-4 mr-2" />
              <SelectValue placeholder="Property" />
            </SelectTrigger>
            <SelectContent>
              {propertyTypes.map((property) => (
                <SelectItem key={property} value={property}>{property}</SelectItem>
              ))}
            </SelectContent>
          </Select>
          
          <Button variant="outline" size="icon">
            <SlidersHorizontal className="h-4 w-4" />
          </Button>
        </div>
      </div>
      
      <Card>
        <CardHeader>
          <CardTitle>Property Explorer</CardTitle>
          <CardDescription>
            Explore molecular properties across the database
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs defaultValue="chart">
            <TabsList className="mb-6">
              <TabsTrigger value="chart">Chart View</TabsTrigger>
              <TabsTrigger value="table">Table View</TabsTrigger>
              <TabsTrigger value="correlation">Correlation Analysis</TabsTrigger>
            </TabsList>
            
            <TabsContent value="chart">
              {selectedProperty ? (
                <div className="space-y-4">
                  <h3 className="text-lg font-medium">{selectedProperty}</h3>
                  
                  {/* Simple visualization placeholder */}
                  <div className="h-[400px] bg-muted/30 rounded-md flex items-center justify-center">
                    <p className="text-muted-foreground">
                      Interactive chart visualization would appear here using Plotly.js
                    </p>
                  </div>
                  
                  <div className="flex justify-end mt-4">
                    <Button variant="outline" className="mr-2">Export Data</Button>
                    <Button>Generate Report</Button>
                  </div>
                </div>
              ) : (
                <div className="text-center py-16">
                  <p className="text-muted-foreground">
                    Select a property from the dropdown above to visualize data.
                  </p>
                </div>
              )}
            </TabsContent>
            
            <TabsContent value="table">
              {selectedProperty ? (
                <div className="space-y-4">
                  <h3 className="text-lg font-medium">{selectedProperty}</h3>
                  
                  {/* Table visualization placeholder */}
                  <div className="border rounded-md overflow-hidden">
                    <table className="w-full">
                      <thead>
                        <tr className="bg-muted">
                          <th className="px-4 py-2 text-left">Molecule</th>
                          <th className="px-4 py-2 text-left">{selectedProperty}</th>
                        </tr>
                      </thead>
                      <tbody>
                        {mockData.labels.map((label, index) => (
                          <tr key={label} className="border-t">
                            <td className="px-4 py-2">{label}</td>
                            <td className="px-4 py-2">{mockData.values[index]} {mockData.unit}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                  
                  <div className="flex justify-end mt-4">
                    <Button variant="outline" className="mr-2">Export Data</Button>
                    <Button>Generate Report</Button>
                  </div>
                </div>
              ) : (
                <div className="text-center py-16">
                  <p className="text-muted-foreground">
                    Select a property from the dropdown above to view data in table format.
                  </p>
                </div>
              )}
            </TabsContent>
            
            <TabsContent value="correlation">
              <div className="space-y-4">
                <div className="flex justify-between items-center">
                  <h3 className="text-lg font-medium">Correlation Analysis</h3>
                  
                  <div className="flex gap-2">
                    <Select>
                      <SelectTrigger className="w-[180px]">
                        <SelectValue placeholder="Property X" />
                      </SelectTrigger>
                      <SelectContent>
                        {propertyTypes.map((property) => (
                          <SelectItem key={property} value={property}>{property}</SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                    
                    <Select>
                      <SelectTrigger className="w-[180px]">
                        <SelectValue placeholder="Property Y" />
                      </SelectTrigger>
                      <SelectContent>
                        {propertyTypes.map((property) => (
                          <SelectItem key={property} value={property}>{property}</SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                  </div>
                </div>
                
                {/* Correlation visualization placeholder */}
                <div className="h-[400px] bg-muted/30 rounded-md flex items-center justify-center">
                  <p className="text-muted-foreground">
                    Correlation scatter plot would appear here when two properties are selected.
                  </p>
                </div>
                
                <div className="flex justify-end mt-4">
                  <Button variant="outline" className="mr-2">Export Data</Button>
                  <Button>Generate Report</Button>
                </div>
              </div>
            </TabsContent>
          </Tabs>
        </CardContent>
      </Card>
    </div>
  )
}