import React from 'react';
import Link from 'next/link';
import { notFound } from 'next/navigation';
import type { Metadata } from 'next';
import ProtocolDetail from '../../../features/protocols/components/protocol-detail';
import { getProtocolById } from '../../../features/protocols/actions/protocol-actions';

interface ProtocolDetailPageProps {
  params: {
    id: string;
  };
}

export async function generateMetadata({ params }: ProtocolDetailPageProps): Promise<Metadata> {
  const protocol = await getProtocolById(params.id);
  
  if (!protocol) {
    return {
      title: 'Protocol Not Found | CryoProtect',
      description: 'The requested protocol could not be found',
    };
  }
  
  return {
    title: `${protocol.name} | CryoProtect`,
    description: protocol.description,
  };
}

export default async function ProtocolDetailPage({ params }: ProtocolDetailPageProps) {
  const protocol = await getProtocolById(params.id);
  
  if (!protocol) {
    notFound();
  }
  
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
      
      <ProtocolDetail protocol={protocol} />
    </div>
  );
}