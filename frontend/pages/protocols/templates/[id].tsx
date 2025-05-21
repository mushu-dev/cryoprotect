/**
 * Edit Protocol Template Page
 */

import React from 'react';
import Head from 'next/head';
import { useRouter } from 'next/router';
import { TemplateEditor } from '../../../src/features/protocols/components/template-editor';
import Layout from '../../../src/components/layout';

export default function EditProtocolTemplatePage() {
  const router = useRouter();
  const { id } = router.query;
  
  return (
    <Layout>
      <Head>
        <title>Edit Protocol Template | CryoProtect</title>
        <meta 
          name="description" 
          content="Edit an existing protocol template for cryopreservation research" 
        />
      </Head>
      
      <main className="container mx-auto px-4 py-8">
        {id && <TemplateEditor templateId={id as string} />}
      </main>
    </Layout>
  );
}