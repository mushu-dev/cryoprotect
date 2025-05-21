import React from 'react';
import Layout from '../components/Layout';

export default function About() {
  return (
    <Layout title="About">
      <div className="page-header">
        <h1 className="title">About CryoProtect</h1>
        <p className="description">
          CryoProtect is a comprehensive platform for cryoprotectant research, analysis, and experimental tracking.
        </p>
      </div>
      
      <div className="card">
        <h2 className="subtitle">Our Mission</h2>
        <p>
          Our mission is to accelerate cryopreservation research and applications by providing integrated tools for molecular analysis, protocol development, and experimental tracking. By centralizing these capabilities, we aim to enable researchers to develop more effective preservation strategies for biological materials.
        </p>
      </div>

      <div className="card">
        <h2 className="subtitle">Key Features</h2>
        <ul className="list">
          <li><strong>Molecular Database:</strong> Comprehensive information on cryoprotectant molecules and their properties</li>
          <li><strong>Mixture Analysis:</strong> Tools to design and analyze cryoprotectant mixtures</li>
          <li><strong>Protocol Management:</strong> Create, optimize, and share preservation protocols</li>
          <li><strong>Experiment Tracking:</strong> Document experiments and analyze results</li>
          <li><strong>Visualization Tools:</strong> Interactive molecular visualization and data presentation</li>
        </ul>
      </div>

      <div className="card">
        <h2 className="subtitle">Research Applications</h2>
        <p>
          CryoProtect supports a wide range of preservation applications, including:
        </p>
        <ul className="list">
          <li>Cell line preservation</li>
          <li>Tissue banking</li>
          <li>Organ preservation</li>
          <li>Protein stabilization</li>
          <li>Seed and plant material storage</li>
          <li>Biopharmaceutical development</li>
        </ul>
      </div>

      <div className="card">
        <h2 className="subtitle">Contact</h2>
        <p>
          For more information about CryoProtect or to provide feedback, please contact us at:
        </p>
        <p><strong>Email:</strong> info@cryoprotect-research.org</p>
        <p><strong>Twitter:</strong> @CryoProtectApp</p>
      </div>
    </Layout>
  );
}