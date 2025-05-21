import React from 'react';
import Head from 'next/head';
import { 
  Clipboard, 
  Check,
  AlertCircle,
  Info,
  ArrowRight,
  Trash,
  Plus
} from 'lucide-react';

import DashboardLayout from '../components/layouts/dashboard-layout';
import { Button } from '../components/ui/button.jsx';
import { Card, CardHeader, CardTitle, CardDescription, CardContent, CardFooter } from '../components/ui/card.jsx';
import { Separator } from '../components/ui/separator.jsx';
import { Alert, AlertDescription, AlertTitle } from '../components/ui/alert.jsx';
import { Badge } from '../components/ui/badge.jsx';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../components/ui/tabs.jsx';
import { Switch } from '../components/ui/switch.jsx';
import { Input } from '../components/ui/input.jsx';
import { Label } from '../components/ui/label.jsx';
import { Textarea } from '../components/ui/textarea.jsx';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '../components/ui/select.jsx';
import { Dialog, DialogTrigger, DialogContent, DialogHeader, DialogTitle, DialogDescription, DialogFooter } from '../components/ui/dialog.jsx';
import { Avatar, AvatarImage, AvatarFallback } from '../components/ui/avatar.jsx';
import { Skeleton } from '../components/ui/skeleton.jsx';

// Component Preview Section
const ComponentSection = ({ title, description, children }) => (
  <section className="mb-12">
    <h2 className="text-2xl font-bold mb-2">{title}</h2>
    <p className="text-muted-foreground mb-6">{description}</p>
    
    <Card>
      <CardContent className="py-6">
        {children}
      </CardContent>
    </Card>
  </section>
);

export default function ComponentsPage() {
  return (
    <DashboardLayout>
      <Head>
        <title>Component Library | CryoProtect</title>
        <meta name="description" content="CryoProtect UI component library and design system" />
      </Head>

      <div className="flex flex-col gap-8">
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Component Library</h1>
          <p className="text-muted-foreground mt-2">
            A documentation of UI components available in the CryoProtect application.
          </p>
          <Separator className="my-6" />
        </div>

        {/* Button Examples */}
        <ComponentSection
          title="Buttons"
          description="Interactive elements for triggering actions."
        >
          <div className="flex flex-wrap gap-4">
            <Button>Default Button</Button>
            <Button variant="secondary">Secondary Button</Button>
            <Button variant="destructive">Destructive Button</Button>
            <Button variant="outline">Outline Button</Button>
            <Button variant="ghost">Ghost Button</Button>
            <Button variant="link">Link Button</Button>
          </div>
          <Separator className="my-4" />
          <div className="flex flex-wrap gap-4">
            <Button size="sm">Small Button</Button>
            <Button>Default Size</Button>
            <Button size="lg">Large Button</Button>
            <Button size="icon"><Plus className="h-4 w-4" /></Button>
          </div>
          <Separator className="my-4" />
          <div className="flex flex-wrap gap-4">
            <Button disabled>Disabled Button</Button>
            <Button variant="secondary" disabled>Disabled Secondary</Button>
            <Button className="flex items-center gap-2">
              Continue <ArrowRight className="h-4 w-4" />
            </Button>
          </div>
        </ComponentSection>

        {/* Card Examples */}
        <ComponentSection
          title="Cards"
          description="Containers for organizing content and data."
        >
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <Card>
              <CardHeader>
                <CardTitle>Simple Card</CardTitle>
                <CardDescription>Basic card with header and content.</CardDescription>
              </CardHeader>
              <CardContent>
                <p>This is a simple card component for displaying content in a contained space.</p>
              </CardContent>
              <CardFooter>
                <Button>Action</Button>
              </CardFooter>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle>Data Card</CardTitle>
                <CardDescription>Card with data representation.</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-2">
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Success Rate</span>
                    <span className="font-medium">85%</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Total Entries</span>
                    <span className="font-medium">42</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Last Updated</span>
                    <span className="font-medium">12 hours ago</span>
                  </div>
                </div>
              </CardContent>
              <CardFooter className="flex justify-between">
                <Button variant="ghost">View Details</Button>
                <Button variant="outline">Export</Button>
              </CardFooter>
            </Card>
          </div>
        </ComponentSection>
        
        {/* Alert Examples */}
        <ComponentSection
          title="Alerts"
          description="Informational messages for user feedback."
        >
          <div className="space-y-4">
            <Alert>
              <Info className="h-4 w-4" />
              <AlertTitle>Information</AlertTitle>
              <AlertDescription>
                Standard information alert with neutral styling.
              </AlertDescription>
            </Alert>
            
            <Alert variant="destructive">
              <AlertCircle className="h-4 w-4" />
              <AlertTitle>Error</AlertTitle>
              <AlertDescription>
                This action cannot be completed due to missing parameters.
              </AlertDescription>
            </Alert>
            
            <Alert className="bg-green-50 text-green-800 border-green-200">
              <Check className="h-4 w-4" />
              <AlertTitle>Success</AlertTitle>
              <AlertDescription>
                Your experiment has been successfully saved.
              </AlertDescription>
            </Alert>
          </div>
        </ComponentSection>

        {/* Badges */}
        <ComponentSection
          title="Badges"
          description="Small labels for status and metadata."
        >
          <div className="flex flex-wrap gap-4">
            <Badge>Default</Badge>
            <Badge variant="secondary">Secondary</Badge>
            <Badge variant="destructive">Destructive</Badge>
            <Badge variant="outline">Outline</Badge>
            <Badge className="bg-blue-500">Custom Blue</Badge>
            <Badge className="bg-amber-500">Custom Amber</Badge>
            <Badge className="bg-green-500">Custom Green</Badge>
            <Badge className="bg-purple-500">Custom Purple</Badge>
          </div>
        </ComponentSection>

        {/* Tabs */}
        <ComponentSection
          title="Tabs"
          description="Organize content into multiple sections."
        >
          <Tabs defaultValue="overview">
            <TabsList className="mb-4">
              <TabsTrigger value="overview">Overview</TabsTrigger>
              <TabsTrigger value="experiments">Experiments</TabsTrigger>
              <TabsTrigger value="properties">Properties</TabsTrigger>
              <TabsTrigger value="analysis">Analysis</TabsTrigger>
            </TabsList>
            <TabsContent value="overview">
              <div className="p-4 border rounded-md">
                <h3 className="text-lg font-medium mb-2">Overview Tab</h3>
                <p>This tab displays the general overview information about the molecule.</p>
              </div>
            </TabsContent>
            <TabsContent value="experiments">
              <div className="p-4 border rounded-md">
                <h3 className="text-lg font-medium mb-2">Experiments Tab</h3>
                <p>View experiments related to this molecule or mixture.</p>
              </div>
            </TabsContent>
            <TabsContent value="properties">
              <div className="p-4 border rounded-md">
                <h3 className="text-lg font-medium mb-2">Properties Tab</h3>
                <p>Physical and chemical properties are displayed here.</p>
              </div>
            </TabsContent>
            <TabsContent value="analysis">
              <div className="p-4 border rounded-md">
                <h3 className="text-lg font-medium mb-2">Analysis Tab</h3>
                <p>Data analysis and insights about performance.</p>
              </div>
            </TabsContent>
          </Tabs>
        </ComponentSection>

        {/* Form Controls */}
        <ComponentSection
          title="Form Controls"
          description="UI elements for building forms and collecting user input."
        >
          <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
            <div className="space-y-6">
              <div className="space-y-2">
                <Label htmlFor="name">Name</Label>
                <Input id="name" placeholder="Enter your name" />
              </div>
              
              <div className="space-y-2">
                <Label htmlFor="email">Email</Label>
                <Input id="email" type="email" placeholder="Enter your email" />
              </div>
              
              <div className="space-y-2">
                <Label htmlFor="category">Category</Label>
                <Select>
                  <SelectTrigger id="category">
                    <SelectValue placeholder="Select a category" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="cryoprotectant">Cryoprotectant</SelectItem>
                    <SelectItem value="penetrating">Penetrating Agent</SelectItem>
                    <SelectItem value="non-penetrating">Non-Penetrating Agent</SelectItem>
                    <SelectItem value="mixture">Mixture</SelectItem>
                  </SelectContent>
                </Select>
              </div>
            </div>
            
            <div className="space-y-6">
              <div className="space-y-2">
                <Label htmlFor="notes">Notes</Label>
                <Textarea id="notes" placeholder="Additional notes or observations..." />
              </div>
              
              <div className="flex items-center gap-2 pt-2">
                <Switch id="notifications" />
                <Label htmlFor="notifications">Enable notifications</Label>
              </div>
              
              <div className="pt-4">
                <Button className="mr-2">Submit</Button>
                <Button variant="outline">Cancel</Button>
              </div>
            </div>
          </div>
        </ComponentSection>

        {/* Dialog */}
        <ComponentSection
          title="Dialogs"
          description="Modal interfaces for focused interactions."
        >
          <div className="flex gap-4">
            <Dialog>
              <DialogTrigger asChild>
                <Button>Open Dialog</Button>
              </DialogTrigger>
              <DialogContent>
                <DialogHeader>
                  <DialogTitle>Confirm Action</DialogTitle>
                  <DialogDescription>
                    This action cannot be undone. Are you sure you want to continue?
                  </DialogDescription>
                </DialogHeader>
                <div className="py-4">
                  <p>Additional information or content can be placed here.</p>
                </div>
                <DialogFooter>
                  <Button variant="outline">Cancel</Button>
                  <Button>Confirm</Button>
                </DialogFooter>
              </DialogContent>
            </Dialog>

            <Dialog>
              <DialogTrigger asChild>
                <Button variant="destructive">Delete Item</Button>
              </DialogTrigger>
              <DialogContent>
                <DialogHeader>
                  <DialogTitle>Delete Item</DialogTitle>
                  <DialogDescription>
                    This will permanently delete the item and all associated data.
                  </DialogDescription>
                </DialogHeader>
                <DialogFooter>
                  <Button variant="outline">Cancel</Button>
                  <Button variant="destructive">
                    <Trash className="mr-2 h-4 w-4" />
                    Delete
                  </Button>
                </DialogFooter>
              </DialogContent>
            </Dialog>
          </div>
        </ComponentSection>

        {/* Avatars and Skeletons */}
        <ComponentSection
          title="Avatars & Loading States"
          description="User representations and loading indicators."
        >
          <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
            <div>
              <h3 className="text-lg font-medium mb-4">Avatars</h3>
              <div className="flex gap-4">
                <Avatar>
                  <AvatarFallback>JD</AvatarFallback>
                </Avatar>
                <Avatar>
                  <AvatarFallback>AB</AvatarFallback>
                </Avatar>
                <Avatar>
                  <AvatarFallback>CP</AvatarFallback>
                </Avatar>
                <Avatar>
                  <AvatarFallback>
                    <User className="h-4 w-4" />
                  </AvatarFallback>
                </Avatar>
              </div>
            </div>
            
            <div>
              <h3 className="text-lg font-medium mb-4">Loading States</h3>
              <div className="space-y-3">
                <div className="flex items-center space-x-4">
                  <Skeleton className="h-12 w-12 rounded-full" />
                  <div className="space-y-2">
                    <Skeleton className="h-4 w-[250px]" />
                    <Skeleton className="h-4 w-[200px]" />
                  </div>
                </div>
                <Skeleton className="h-4 w-full" />
                <Skeleton className="h-4 w-full" />
                <Skeleton className="h-4 w-3/4" />
              </div>
            </div>
          </div>
        </ComponentSection>
      </div>
    </DashboardLayout>
  );
}