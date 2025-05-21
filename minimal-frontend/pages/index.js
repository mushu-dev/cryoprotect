import React from 'react';
import Layout from '../components/Layout';

export default function Home() {
  return (
    <Layout title="Home">
      <h1 className="title">Welcome to CryoProtect</h1>
      
      <p className="description">
        CryoProtect is a platform for cryoprotectant research, analysis, and experimental tracking. 
        We're currently working on enhancements to our platform to better serve the scientific community.
      </p>

      <div className="card">
        <h2 className="subtitle">Coming Soon</h2>
        <ul className="list">
          <li>Enhanced molecular visualization</li>
          <li>Improved experimental data tracking</li>
          <li>Protocol and experiment sharing</li>
          <li>Advanced analysis tools</li>
          <li>Integration with laboratory equipment</li>
        </ul>
      </div>
      
      <div className="card">
        <h2 className="subtitle">Available Features</h2>
        <div className="feature-grid">
          <div className="feature-card">
            <h3>Dashboard</h3>
            <p>View statistics and an overview of your data.</p>
            <a href="/dashboard" className="feature-link">Go to Dashboard</a>
          </div>
          
          <div className="feature-card">
            <h3>Convex Integration</h3>
            <p>Explore our real-time database integration with Convex.</p>
            <a href="/convex" className="feature-link">View Convex Demo</a>
          </div>
          
          <div className="feature-card">
            <h3>Molecules</h3>
            <p>Browse available molecules in our database.</p>
            <a href="/molecules" className="feature-link">View Molecules</a>
          </div>
          
          <div className="feature-card">
            <h3>Mixtures</h3>
            <p>Explore cryoprotectant mixtures and their compositions.</p>
            <a href="/mixtures" className="feature-link">View Mixtures</a>
          </div>
          
          <div className="feature-card">
            <h3>Backend Connections</h3>
            <p>View the status of connections to our backend services.</p>
            <a href="/connections" className="feature-link">Check Connections</a>
          </div>
        </div>
      </div>
      
      <style jsx>{`
        .feature-grid {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
          gap: 20px;
          margin-top: 15px;
        }
        
        .feature-card {
          padding: 15px;
          border: 1px solid #eee;
          border-radius: 5px;
          background-color: #f9f9f9;
        }
        
        .feature-card h3 {
          margin-top: 0;
          color: #3182ce;
        }
        
        .feature-link {
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
        
        .feature-link:hover {
          background-color: #2c5282;
        }
      `}</style>
    </Layout>
  );
}