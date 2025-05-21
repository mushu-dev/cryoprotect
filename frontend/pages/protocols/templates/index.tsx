/**
 * Protocol Templates List Page
 */

import React from 'react';
import Head from 'next/head';
import { TemplateManagement } from '../../../src/features/protocols/components/template-management';
import Layout from '../../../src/components/layout';

export default function ProtocolTemplatesPage() {
  return (
    <Layout>
      <Head>
        <title>Protocol Templates | CryoProtect</title>
        <meta 
          name="description" 
          content="Manage and use standardized protocol templates for cryopreservation research" 
        />
      </Head>
      
      <main className="container mx-auto px-4 py-8">
        <TemplateManagement />
      </main>
    </Layout>
  );
}