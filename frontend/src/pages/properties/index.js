import React from 'react';
import Head from 'next/head';

export default function PropertiesPage() {
  return (
    <>
      <Head>
        <title>Properties - CryoProtect</title>
        <meta name="description" content="Explore cryoprotectant properties and measurements" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Properties</h1>
        
        <div className="prose max-w-none">
          <p>
            This page provides access to the physical and chemical properties of cryoprotectant molecules 
            in our database. You can explore different property categories and their relationships.
          </p>
          
          <h2 className="text-2xl font-semibold mt-8 mb-4">Property Categories</h2>
          
          <ul className="list-disc pl-5 space-y-2">
            <li>Physical Properties (melting point, boiling point, density)</li>
            <li>Chemical Properties (molecular weight, LogP, pKa)</li>
            <li>Cryoprotective Properties (glass transition temperature, vitrification)</li>
            <li>Toxicity Indicators (cytotoxicity, genotoxicity)</li>
            <li>Experimental Measurements (cell viability, membrane permeability)</li>
          </ul>
          
          <div className="bg-blue-50 border-l-4 border-blue-500 p-4 my-6">
            <p className="text-blue-700">
              Detailed property exploration features are coming soon. 
              Check back for interactive property visualization and comparison tools.
            </p>
          </div>
        </div>
      </div>
    </>
  );
}