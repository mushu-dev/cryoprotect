/**
 * Create Protocol Template Page
 */

import React from 'react';
import Head from 'next/head';
import { TemplateEditor } from '../../../src/features/protocols/components/template-editor';
import Layout from '../../../src/components/layout';

export default function CreateProtocolTemplatePage() {
  return (
    <Layout>
      <Head>
        <title>Create Protocol Template | CryoProtect</title>
        <meta 
          name="description" 
          content="Create a new protocol template for cryopreservation research" 
        />
      </Head>
      
      <main className="container mx-auto px-4 py-8">
        <TemplateEditor />
      </main>
    </Layout>
  );
}