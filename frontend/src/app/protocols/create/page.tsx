import React from 'react';
import { Metadata } from 'next';
import { Button } from '../../../components/ui/button';
import { ProtocolBuilder } from '../../../features/protocols/components/protocol-builder';
import Link from 'next/link';

export const metadata: Metadata = {
  title: 'Create Protocol | CryoProtect',
  description: 'Create a new standardized protocol for cryopreservation experiments',
};

export default function CreateProtocolPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Create New Protocol</h1>
          <p className="text-muted-foreground">
            Design a standardized protocol for cryopreservation procedures
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Link href="/protocols" passHref>
            <Button variant="outline">Cancel</Button>
          </Link>
        </div>
      </div>

      <ProtocolBuilderWrapper />
    </div>
  );
}

// Client component wrapper
function ProtocolBuilderWrapper() {
  // In a real implementation, this would handle form submission logic
  const handleSave = (protocolData: any) => {
    console.log('Saving protocol:', protocolData);
    // Would redirect to the new protocol page or protocols list
  };

  const handleCancel = () => {
    // Would navigate back to protocols list
  };

  return (
    <div className="bg-card rounded-lg shadow-sm overflow-hidden">
      <ProtocolBuilder 
        onSave={handleSave}
        onCancel={handleCancel}
      />
    </div>
  );
}