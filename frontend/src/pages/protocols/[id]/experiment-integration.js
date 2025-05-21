import React from 'react';
import { useRouter } from 'next/router';
import Link from 'next/link';
import Head from 'next/head';
import { useConvex } from 'convex/react';
import { ProtocolExperimentIntegration } from '../../../features/protocols/components';
import { ChevronLeft, ExperimentIcon, FileText, Home, ScienceIcon } from 'lucide-react';

export default function ProtocolExperimentIntegrationPage() {
  const router = useRouter();
  const { id } = router.query;
  const convex = useConvex();
  
  if (router.isFallback || !id) {
    return (
      <div className="container mx-auto px-4 py-8">
        <div className="flex justify-center items-center h-64">
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-primary mx-auto mb-4"></div>
            <p className="text-muted-foreground">Loading protocol integration...</p>
          </div>
        </div>
      </div>
    );
  }
  
  return (
    <>
      <Head>
        <title>Protocol-Experiment Integration | CryoProtect</title>
        <meta name="description" content="Integrate protocols with experiments for standardized workflows" />
      </Head>
      
      <div className="container mx-auto px-4 py-8">
        {/* Breadcrumb Navigation */}
        <nav className="flex mb-6" aria-label="Breadcrumb">
          <ol className="inline-flex items-center space-x-1 md:space-x-3">
            <li className="inline-flex items-center">
              <Link href="/" className="inline-flex items-center text-sm font-medium text-muted-foreground hover:text-foreground">
                <Home className="w-4 h-4 mr-2.5" />
                Home
              </Link>
            </li>
            <li>
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <Link href="/protocols" className="ml-1 text-sm font-medium text-muted-foreground hover:text-foreground md:ml-2">
                  <span className="flex items-center">
                    <FileText className="w-4 h-4 mr-1" />
                    Protocols
                  </span>
                </Link>
              </div>
            </li>
            <li>
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <Link href={`/protocols/${id}`} className="ml-1 text-sm font-medium text-muted-foreground hover:text-foreground md:ml-2">
                  Protocol Details
                </Link>
              </div>
            </li>
            <li aria-current="page">
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <span className="ml-1 text-sm font-medium text-foreground md:ml-2 flex items-center">
                  <ScienceIcon className="w-4 h-4 mr-1" />
                  Experiments
                </span>
              </div>
            </li>
          </ol>
        </nav>
        
        {/* Header */}
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-6">
          <h1 className="text-3xl font-bold mb-4 md:mb-0">Protocol-Experiment Integration</h1>
          <div className="flex space-x-2">
            <Link href={`/protocols/${id}`} className="inline-flex items-center justify-center rounded-md bg-muted px-4 py-2 text-sm font-medium text-muted-foreground shadow-sm hover:bg-muted/80">
              <ChevronLeft className="mr-2 h-4 w-4" />
              Back to Protocol
            </Link>
            <Link href="/experiments" className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow-sm hover:bg-secondary/80">
              <ScienceIcon className="mr-2 h-4 w-4" />
              All Experiments
            </Link>
          </div>
        </div>
        
        {/* Main content */}
        <ProtocolExperimentIntegration 
          protocolId={id}
          convexClient={convex}
        />
      </div>
    </>
  );
}