import React from 'react';
import Layout from '../components/Layout';

export default function Protocols() {
  return (
    <Layout title="Protocols">
      <div className="page-header">
        <h1 className="title">Protocols</h1>
        <p className="description">
          Create, manage, and share cryopreservation protocols. Step-by-step procedures, equipment configurations, and optimization tools will be available soon.
        </p>
      </div>
      
      <div className="card">
        <h2 className="subtitle">Coming Features</h2>
        <ul className="list">
          <li>Protocol Step Editor</li>
          <li>Equipment and Parameter Management</li>
          <li>Version Control and Comparison</li>
          <li>Protocol Validation Tools</li>
          <li>Template Library</li>
        </ul>
      </div>

      <div className="card">
        <h2 className="subtitle">Protocol Types</h2>
        <ul className="list">
          <li>Cell Cryopreservation</li>
          <li>Tissue Preservation</li>
          <li>Protein Stabilization</li>
          <li>Organoid Preservation</li>
          <li>Custom Protocols</li>
        </ul>
      </div>
    </Layout>
  );
}