import React from 'react';
import Layout from '../components/Layout';
import ConnectionStatus from '../src/components/convex/ConnectionStatus';

export default function ConnectionsPage() {
  return (
    <Layout title="Backend Connections">
      <div className="page-header">
        <h1 className="title">Backend Connections</h1>
        <p className="description">
          This page shows the status of connections to our various backend services.
        </p>
      </div>

      <div className="card">
        <ConnectionStatus />
      </div>

      <div className="card" style={{ marginTop: '30px' }}>
        <h2 className="subtitle">Backend Services</h2>
        <div className="backend-services">
          <div className="service-item">
            <h3>Heroku API</h3>
            <p>Main Flask API that serves molecules, mixtures, and experimental data.</p>
            <p><strong>URL:</strong> https://cryoprotect-8030e4025428.herokuapp.com/api</p>
            <p><strong>Endpoints:</strong></p>
            <ul>
              <li>/molecules - Get all molecules</li>
              <li>/molecules/:id - Get molecule by ID</li>
              <li>/mixtures - Get all mixtures</li>
              <li>/mixtures/:id - Get mixture by ID</li>
            </ul>
          </div>

          <div className="service-item">
            <h3>Fly.io RDKit Service</h3>
            <p>Specialized service for molecular calculations and visualizations.</p>
            <p><strong>URL:</strong> https://cryoprotect-rdkit.fly.dev</p>
            <p><strong>Endpoints:</strong></p>
            <ul>
              <li>/calculate - Calculate molecular properties</li>
              <li>/depict - Generate 2D molecular depictions</li>
              <li>/similarity - Calculate molecular similarity</li>
              <li>/health - Health check endpoint</li>
            </ul>
          </div>

          <div className="service-item">
            <h3>Convex Database</h3>
            <p>Real-time database with automatic synchronization.</p>
            <p><strong>URL:</strong> https://primary-meerkat-478.convex.cloud</p>
            <p><strong>Features:</strong></p>
            <ul>
              <li>Real-time data updates</li>
              <li>Full TypeScript integration</li>
              <li>Automatic client-side caching</li>
              <li>Authentication integration</li>
            </ul>
            <p><a href="/convex" className="demo-link">View Convex Demo Page</a></p>
          </div>
        </div>
      </div>
      
      <style jsx>{`
        .service-item {
          margin-bottom: 20px;
          padding: 15px;
          border: 1px solid #eee;
          border-radius: 5px;
        }
        
        .backend-services {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
          gap: 20px;
        }
        
        .demo-link {
          display: inline-block;
          margin-top: 10px;
          padding: 8px 16px;
          background-color: #3182ce;
          color: white;
          border-radius: 4px;
          text-decoration: none;
          font-weight: 500;
          transition: background-color 0.2s ease;
        }
        
        .demo-link:hover {
          background-color: #2c5282;
        }
      `}</style>
    </Layout>
  );
}