import React from 'react';
import { Metadata } from 'next';
import Link from 'next/link';
import { Button } from '../../components/ui/button';
import { Card } from '../../components/ui/card';
import { ProtocolCard } from '../../features/protocols/components/protocol-card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../components/ui/tabs';
import { Input } from '../../components/ui/input';

export const metadata: Metadata = {
  title: 'Protocols | CryoProtect',
  description: 'Browse and manage standardized protocols for cryopreservation experiments',
};

// This would be a server component in a real implementation
export default function ProtocolsPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Protocols</h1>
          <p className="text-muted-foreground">
            Standardized protocols for reproducible cryopreservation procedures
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Link href="/protocols/create" passHref>
            <Button>Create New Protocol</Button>
          </Link>
        </div>
      </div>

      <ProtocolsContent />
    </div>
  );
}

// Client component for dynamic content
function ProtocolsContent() {
  // This is a client component placeholder
  // In a real implementation, this would use the useProtocols hook
  
  // Mock data for display
  const myProtocols = [
    {
      id: '1',
      name: 'Standard Cell Freezing',
      description: 'Basic protocol for freezing mammalian cells using DMSO',
      version: '1.2.0',
      created_at: '2023-09-01',
      updated_at: '2023-10-15',
      step_count: 8,
      duration: '2 hours',
      is_template: false,
      author: { name: 'John Doe', id: 'user1' }
    },
    {
      id: '2',
      name: 'Vitrification Protocol',
      description: 'Advanced vitrification method with controlled cooling rates',
      version: '2.0.1',
      created_at: '2023-08-12',
      updated_at: '2023-10-18',
      step_count: 12,
      duration: '3.5 hours',
      is_template: false,
      author: { name: 'John Doe', id: 'user1' }
    }
  ];
  
  const templateProtocols = [
    {
      id: '3',
      name: 'Standard Slow Freezing',
      description: 'Industry standard slow freezing protocol for various cell types',
      version: '3.1.0',
      created_at: '2023-06-15',
      updated_at: '2023-09-20',
      step_count: 10,
      duration: '2.5 hours',
      is_template: true,
      author: { name: 'System', id: 'system' }
    },
    {
      id: '4',
      name: 'Rapid Cooling Method',
      description: 'Quick cooling protocol for small sample volumes',
      version: '1.5.2',
      created_at: '2023-07-22',
      updated_at: '2023-08-30',
      step_count: 6,
      duration: '1 hour',
      is_template: true,
      author: { name: 'System', id: 'system' }
    },
    {
      id: '5',
      name: 'Plant Tissue Cryopreservation',
      description: 'Specialized protocol for preserving plant tissue samples',
      version: '2.2.0',
      created_at: '2023-05-18',
      updated_at: '2023-09-10',
      step_count: 14,
      duration: '4 hours',
      is_template: true,
      author: { name: 'Jane Smith', id: 'user2' }
    }
  ];

  return (
    <>
      <div className="mb-8">
        <div className="max-w-md mx-auto mb-8">
          <Input 
            type="search" 
            placeholder="Search protocols..." 
            className="w-full" 
          />
        </div>

        <Tabs defaultValue="all">
          <TabsList className="mb-6">
            <TabsTrigger value="all">All Protocols</TabsTrigger>
            <TabsTrigger value="my">My Protocols</TabsTrigger>
            <TabsTrigger value="templates">Templates</TabsTrigger>
          </TabsList>
          
          <TabsContent value="all">
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {[...myProtocols, ...templateProtocols].map(protocol => (
                <Link key={protocol.id} href={`/protocols/${protocol.id}`} passHref>
                  <ProtocolCard 
                    protocol={protocol} 
                    onUseTemplate={protocol.is_template ? () => {} : undefined}
                  />
                </Link>
              ))}
            </div>
          </TabsContent>
          
          <TabsContent value="my">
            {myProtocols.length > 0 ? (
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
                {myProtocols.map(protocol => (
                  <Link key={protocol.id} href={`/protocols/${protocol.id}`} passHref>
                    <ProtocolCard protocol={protocol} />
                  </Link>
                ))}
              </div>
            ) : (
              <div className="text-center py-12 bg-card rounded-lg">
                <p className="mb-4">You haven't created any protocols yet.</p>
                <Link href="/protocols/create" passHref>
                  <Button>Create Your First Protocol</Button>
                </Link>
              </div>
            )}
          </TabsContent>
          
          <TabsContent value="templates">
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {templateProtocols.map(protocol => (
                <ProtocolCard 
                  key={protocol.id} 
                  protocol={protocol} 
                  onUseTemplate={() => {}}
                />
              ))}
            </div>
          </TabsContent>
        </Tabs>
      </div>
    </>
  );
}