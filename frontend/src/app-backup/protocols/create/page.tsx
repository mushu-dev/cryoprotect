import React from 'react';
import type { Metadata } from 'next';
import Link from 'next/link';
import { ProtocolBuilder } from '../../../features/protocols/components/protocol-builder';

export const metadata: Metadata = {
  title: 'Create Protocol | CryoProtect',
  description: 'Create a new standardized protocol for cryopreservation experiments',
};

export default function CreateProtocolPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="mb-6">
        <Link 
          href="/protocols"
          className="text-primary hover:underline flex items-center"
        >
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
            <path d="M19 12H5"/>
            <path d="M12 19l-7-7 7-7"/>
          </svg>
          Back to Protocols
        </Link>
      </div>
      
      <div className="mb-6">
        <h1 className="text-3xl font-bold mb-2">Create New Protocol</h1>
        <p className="text-muted-foreground">
          Design a new standardized protocol with detailed steps, parameters, and materials
        </p>
      </div>
      
      <ProtocolBuilder />
    </div>
  );
}