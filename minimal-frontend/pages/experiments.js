import React from 'react';
import Layout from '../components/Layout';

export default function Experiments() {
  return (
    <Layout title="Experiments">
      <div className="page-header">
        <h1 className="title">Experiments</h1>
        <p className="description">
          Track, analyze, and share cryopreservation experiments. Data collection, visualization, and analysis tools will be available soon.
        </p>
      </div>
      
      <div className="card">
        <h2 className="subtitle">Coming Features</h2>
        <ul className="list">
          <li>Experiment Data Collection</li>
          <li>Results Visualization</li>
          <li>Protocol Execution Tracking</li>
          <li>Comparative Analysis</li>
          <li>Export and Sharing Capabilities</li>
        </ul>
      </div>

      <div className="card">
        <h2 className="subtitle">Data Visualization</h2>
        <p className="description">
          Our platform will offer powerful visualization tools for experimental data, enabling researchers to:
        </p>
        <ul className="list">
          <li>Compare viability across different conditions</li>
          <li>Analyze temperature profiles during freezing and thawing</li>
          <li>Visualize cell recovery rates</li>
          <li>Track experimental reproducibility</li>
          <li>Generate publication-ready figures</li>
        </ul>
      </div>
    </Layout>
  );
}