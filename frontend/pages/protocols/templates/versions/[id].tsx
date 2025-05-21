/**
 * Protocol Template Version History Page
 */

import React from 'react';
import Head from 'next/head';
import { useRouter } from 'next/router';
import { TemplateComparison } from '../../../../src/features/protocols/components/template-comparison';
import Layout from '../../../../src/components/layout';

export default function TemplateVersionHistoryPage() {
  const router = useRouter();
  const { id } = router.query;
  
  return (
    <Layout>
      <Head>
        <title>Template Version History | CryoProtect</title>
        <meta 
          name="description" 
          content="Compare different versions of a protocol template" 
        />
      </Head>
      
      <main className="container mx-auto px-4 py-8">
        {id && <TemplateComparison templateId={id as string} />}
      </main>
    </Layout>
  );
}