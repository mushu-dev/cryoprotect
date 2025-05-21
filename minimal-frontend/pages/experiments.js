import React, { useState, useEffect } from 'react';
import dynamic from 'next/dynamic';
import Layout from '../components/Layout';
import Loading from '../components/Loading';

// Check if Convex is enabled
const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

// Dynamically import Convex components to avoid build issues
const ConvexExperimentsView = dynamic(
  () => import('../src/components/ConvexExperimentsView'),
  { 
    ssr: false,
    loading: () => <Loading />
  }
);

// Mock experiments data for demo
const mockExperiments = [
  {
    id: 1,
    title: 'DMSO Viability Study',
    description: 'Testing cell viability with various DMSO concentrations',
    status: 'completed',
    tissue_type: 'HeLa cells',
    cryoprotectant: 'DMSO',
    concentration: '10%',
    date_created: '2024-01-15',
    date_performed: '2024-01-18',
    results: {
      viability: 85.2,
      recovery_rate: 78.5,
      functionality_score: 82.1
    },
    protocol: {
      cooling_rate: '-1°C/min',
      storage_temp: '-196°C',
      thaw_rate: 'rapid'
    }
  },
  {
    id: 2,
    title: 'Glycerol vs DMSO Comparison',
    description: 'Comparative study of glycerol and DMSO for stem cell preservation',
    status: 'in_progress',
    tissue_type: 'Mesenchymal stem cells',
    cryoprotectant: 'Glycerol + DMSO',
    concentration: '5% + 5%',
    date_created: '2024-01-20',
    date_performed: null,
    results: null,
    protocol: {
      cooling_rate: '-0.5°C/min',
      storage_temp: '-80°C',
      thaw_rate: 'controlled'
    }
  },
  {
    id: 3,
    title: 'Slow Freezing Protocol Optimization',
    description: 'Optimizing cooling rates for improved cell survival',
    status: 'completed',
    tissue_type: 'Primary hepatocytes',
    cryoprotectant: 'Ethylene glycol',
    concentration: '15%',
    date_created: '2024-01-10',
    date_performed: '2024-01-12',
    results: {
      viability: 92.7,
      recovery_rate: 88.3,
      functionality_score: 90.1
    },
    protocol: {
      cooling_rate: '-0.3°C/min',
      storage_temp: '-196°C',
      thaw_rate: 'controlled'
    }
  },
  {
    id: 4,
    title: 'Propylene Glycol Toxicity Assessment',
    description: 'Evaluating toxicity levels of propylene glycol at various concentrations',
    status: 'planned',
    tissue_type: 'Cardiac myocytes',
    cryoprotectant: 'Propylene glycol',
    concentration: 'Variable (5-20%)',
    date_created: '2024-01-25',
    date_performed: null,
    results: null,
    protocol: {
      cooling_rate: '-1°C/min',
      storage_temp: '-80°C',
      thaw_rate: 'rapid'
    }
  }
];

// Standard experiments view component
function StandardExperimentsView({ experiments }) {
  const [searchTerm, setSearchTerm] = useState('');
  const [statusFilter, setStatusFilter] = useState('all');
  const [selectedExperiment, setSelectedExperiment] = useState(null);
  
  // Filter experiments based on search term and status
  const filteredExperiments = React.useMemo(() => {
    if (!experiments) return [];
    
    let filtered = experiments;
    
    // Filter by search term
    if (searchTerm) {
      filtered = filtered.filter(exp => 
        exp.title.toLowerCase().includes(searchTerm.toLowerCase()) ||
        exp.description.toLowerCase().includes(searchTerm.toLowerCase()) ||
        exp.tissue_type.toLowerCase().includes(searchTerm.toLowerCase()) ||
        exp.cryoprotectant.toLowerCase().includes(searchTerm.toLowerCase())
      );
    }
    
    // Filter by status
    if (statusFilter !== 'all') {
      filtered = filtered.filter(exp => exp.status === statusFilter);
    }
    
    return filtered;
  }, [experiments, searchTerm, statusFilter]);
  
  const getStatusColor = (status) => {
    switch (status) {
      case 'completed': return 'status-completed';
      case 'in_progress': return 'status-in-progress';
      case 'planned': return 'status-planned';
      default: return 'status-default';
    }
  };
  
  const getViabilityColor = (viability) => {
    if (!viability) return '';
    if (viability >= 90) return 'excellent';
    if (viability >= 80) return 'good';
    if (viability >= 70) return 'fair';
    return 'poor';
  };

  return (
    <div className="experiments-container">
      <div className="experiments-header">
        <div className="header-content">
          <h1 className="page-title">Experiments Database</h1>
          <p className="page-description">
            Track, analyze, and share cryopreservation experiments with comprehensive data collection and visualization
          </p>
        </div>
        
        {/* Search and filters */}
        <div className="search-section">
          <input
            type="text"
            placeholder="Search experiments by title, description, tissue type, or cryoprotectant..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            className="search-input"
          />
          
          <select
            value={statusFilter}
            onChange={(e) => setStatusFilter(e.target.value)}
            className="filter-select"
          >
            <option value="all">All Status</option>
            <option value="completed">Completed</option>
            <option value="in_progress">In Progress</option>
            <option value="planned">Planned</option>
          </select>
        </div>
        
        {/* Status info */}
        <div className="status-info">
          <span className="status-item">🔬 Experimental Database</span>
          <span className="status-item">📊 Data Analysis</span>
          <span className="status-item">📈 {filteredExperiments?.length || 0} experiments</span>
        </div>
      </div>
      
      <div className="experiments-grid">
        {filteredExperiments && filteredExperiments.length > 0 ? (
          filteredExperiments.map((experiment) => (
            <div key={experiment.id} className="experiment-card">
              <div className="experiment-header">
                <h3 className="experiment-title">{experiment.title}</h3>
                <span className={`experiment-status ${getStatusColor(experiment.status)}`}>
                  {experiment.status.replace('_', ' ')}
                </span>
              </div>
              
              <div className="experiment-details">
                <p className="experiment-description">{experiment.description}</p>
                
                <div className="experiment-metadata">
                  <div className="metadata-row">
                    <span className="metadata-label">Tissue Type:</span>
                    <span className="metadata-value">{experiment.tissue_type}</span>
                  </div>
                  
                  <div className="metadata-row">
                    <span className="metadata-label">Cryoprotectant:</span>
                    <span className="metadata-value">{experiment.cryoprotectant}</span>
                  </div>
                  
                  <div className="metadata-row">
                    <span className="metadata-label">Concentration:</span>
                    <span className="metadata-value">{experiment.concentration}</span>
                  </div>
                  
                  {experiment.date_performed && (
                    <div className="metadata-row">
                      <span className="metadata-label">Date Performed:</span>
                      <span className="metadata-value">{experiment.date_performed}</span>
                    </div>
                  )}
                </div>
              </div>
              
              {/* Protocol section */}
              <div className="protocol-section">
                <h4 className="section-title">Protocol</h4>
                <div className="protocol-details">
                  <div className="protocol-item">
                    <span className="protocol-label">Cooling Rate:</span>
                    <span className="protocol-value">{experiment.protocol.cooling_rate}</span>
                  </div>
                  <div className="protocol-item">
                    <span className="protocol-label">Storage Temp:</span>
                    <span className="protocol-value">{experiment.protocol.storage_temp}</span>
                  </div>
                  <div className="protocol-item">
                    <span className="protocol-label">Thaw Rate:</span>
                    <span className="protocol-value">{experiment.protocol.thaw_rate}</span>
                  </div>
                </div>
              </div>
              
              {/* Results section */}
              {experiment.results ? (
                <div className="results-section">
                  <h4 className="section-title">Results</h4>
                  <div className="results-grid">
                    <div className="result-item">
                      <span className="result-label">Viability:</span>
                      <span className={`result-value ${getViabilityColor(experiment.results.viability)}`}>
                        {experiment.results.viability}%
                      </span>
                    </div>
                    <div className="result-item">
                      <span className="result-label">Recovery:</span>
                      <span className={`result-value ${getViabilityColor(experiment.results.recovery_rate)}`}>
                        {experiment.results.recovery_rate}%
                      </span>
                    </div>
                    <div className="result-item">
                      <span className="result-label">Functionality:</span>
                      <span className={`result-value ${getViabilityColor(experiment.results.functionality_score)}`}>
                        {experiment.results.functionality_score}%
                      </span>
                    </div>
                  </div>
                </div>
              ) : (
                <div className="no-results">
                  <span className="no-results-text">
                    {experiment.status === 'planned' ? 'Experiment not yet performed' : 'Results pending'}
                  </span>
                </div>
              )}
              
              {/* Action buttons */}
              <div className="experiment-actions">
                <button 
                  className="action-btn primary"
                  onClick={() => setSelectedExperiment(experiment)}
                >
                  View Details
                </button>
                {experiment.status === 'completed' && (
                  <button className="action-btn secondary">
                    Download Data
                  </button>
                )}
              </div>
            </div>
          ))
        ) : (
          <div className="no-experiments">
            <p>No experiments found.</p>
            {searchTerm && (
              <p>Try adjusting your search term or filters.</p>
            )}
          </div>
        )}
      </div>
      
      {/* Statistics section */}
      <div className="stats-section">
        <h2 className="stats-title">Experiment Statistics</h2>
        <div className="stats-grid">
          <div className="stat-card">
            <div className="stat-number">{experiments.filter(e => e.status === 'completed').length}</div>
            <div className="stat-label">Completed</div>
          </div>
          <div className="stat-card">
            <div className="stat-number">{experiments.filter(e => e.status === 'in_progress').length}</div>
            <div className="stat-label">In Progress</div>
          </div>
          <div className="stat-card">
            <div className="stat-number">{experiments.filter(e => e.status === 'planned').length}</div>
            <div className="stat-label">Planned</div>
          </div>
          <div className="stat-card">
            <div className="stat-number">
              {experiments.filter(e => e.results && e.results.viability >= 80).length}
            </div>
            <div className="stat-label">High Viability</div>
          </div>
        </div>
      </div>
      
      <style jsx>{`
        .experiments-container {
          max-width: 1200px;
          margin: 0 auto;
          padding: 20px;
        }
        
        .experiments-header {
          margin-bottom: 30px;
        }
        
        .header-content {
          text-align: center;
          margin-bottom: 20px;
        }
        
        .page-title {
          font-size: 2.5rem;
          font-weight: bold;
          color: #1a202c;
          margin-bottom: 10px;
        }
        
        .page-description {
          font-size: 1.1rem;
          color: #4a5568;
          max-width: 700px;
          margin: 0 auto;
        }
        
        .search-section {
          display: flex;
          gap: 15px;
          margin-bottom: 20px;
          align-items: center;
          justify-content: center;
          flex-wrap: wrap;
        }
        
        .search-input {
          flex: 1;
          min-width: 300px;
          max-width: 500px;
          padding: 12px 16px;
          border: 2px solid #e2e8f0;
          border-radius: 8px;
          font-size: 16px;
          outline: none;
          transition: border-color 0.2s;
        }
        
        .search-input:focus {
          border-color: #3182ce;
        }
        
        .filter-select {
          padding: 12px 16px;
          border: 2px solid #e2e8f0;
          border-radius: 8px;
          font-size: 16px;
          outline: none;
          background: white;
          cursor: pointer;
        }
        
        .status-info {
          display: flex;
          justify-content: center;
          gap: 20px;
          font-size: 14px;
          color: #4a5568;
          flex-wrap: wrap;
        }
        
        .experiments-grid {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
          gap: 24px;
          margin-bottom: 40px;
        }
        
        .experiment-card {
          background: white;
          border: 1px solid #e2e8f0;
          border-radius: 12px;
          padding: 24px;
          box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
          transition: transform 0.2s, box-shadow 0.2s;
        }
        
        .experiment-card:hover {
          transform: translateY(-2px);
          box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }
        
        .experiment-header {
          display: flex;
          justify-content: space-between;
          align-items: flex-start;
          margin-bottom: 15px;
          gap: 10px;
        }
        
        .experiment-title {
          font-size: 1.25rem;
          font-weight: 600;
          color: #1a202c;
          margin: 0;
          flex: 1;
        }
        
        .experiment-status {
          padding: 4px 8px;
          border-radius: 4px;
          font-size: 12px;
          font-weight: 500;
          text-transform: capitalize;
          white-space: nowrap;
        }
        
        .status-completed { background: #48bb78; color: white; }
        .status-in-progress { background: #ed8936; color: white; }
        .status-planned { background: #667eea; color: white; }
        .status-default { background: #a0aec0; color: white; }
        
        .experiment-details {
          margin-bottom: 20px;
        }
        
        .experiment-description {
          color: #4a5568;
          margin-bottom: 15px;
          line-height: 1.5;
        }
        
        .experiment-metadata {
          display: grid;
          gap: 8px;
        }
        
        .metadata-row {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 8px;
        }
        
        .metadata-label {
          color: #718096;
          font-size: 14px;
        }
        
        .metadata-value {
          color: #1a202c;
          font-size: 14px;
          font-weight: 500;
        }
        
        .protocol-section, .results-section {
          border-top: 1px solid #e2e8f0;
          padding-top: 16px;
          margin-bottom: 16px;
        }
        
        .section-title {
          font-size: 1rem;
          font-weight: 600;
          color: #1a202c;
          margin: 0 0 12px 0;
        }
        
        .protocol-details {
          display: grid;
          gap: 6px;
        }
        
        .protocol-item {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 8px;
        }
        
        .protocol-label {
          color: #718096;
          font-size: 13px;
        }
        
        .protocol-value {
          color: #1a202c;
          font-size: 13px;
          font-weight: 500;
        }
        
        .results-grid {
          display: grid;
          gap: 8px;
        }
        
        .result-item {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 8px;
        }
        
        .result-label {
          color: #718096;
          font-size: 14px;
        }
        
        .result-value {
          font-size: 14px;
          font-weight: 600;
        }
        
        .result-value.excellent { color: #38a169; }
        .result-value.good { color: #68d391; }
        .result-value.fair { color: #d69e2e; }
        .result-value.poor { color: #e53e3e; }
        
        .no-results {
          border-top: 1px solid #e2e8f0;
          padding-top: 16px;
          margin-bottom: 16px;
        }
        
        .no-results-text {
          color: #718096;
          font-style: italic;
          font-size: 14px;
        }
        
        .experiment-actions {
          display: flex;
          gap: 10px;
          border-top: 1px solid #e2e8f0;
          padding-top: 16px;
        }
        
        .action-btn {
          padding: 8px 16px;
          border-radius: 6px;
          font-size: 14px;
          font-weight: 500;
          cursor: pointer;
          transition: background-color 0.2s;
          border: none;
        }
        
        .action-btn.primary {
          background: #3182ce;
          color: white;
        }
        
        .action-btn.primary:hover {
          background: #2c5282;
        }
        
        .action-btn.secondary {
          background: #f7fafc;
          color: #4a5568;
          border: 1px solid #e2e8f0;
        }
        
        .action-btn.secondary:hover {
          background: #edf2f7;
        }
        
        .no-experiments {
          grid-column: 1 / -1;
          text-align: center;
          color: #4a5568;
          padding: 40px;
        }
        
        .stats-section {
          background: white;
          border: 1px solid #e2e8f0;
          border-radius: 12px;
          padding: 24px;
          margin-top: 30px;
        }
        
        .stats-title {
          font-size: 1.5rem;
          font-weight: 600;
          color: #1a202c;
          margin-bottom: 20px;
          text-align: center;
        }
        
        .stats-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
          gap: 20px;
        }
        
        .stat-card {
          text-align: center;
          padding: 20px;
          background: #f7fafc;
          border-radius: 8px;
        }
        
        .stat-number {
          font-size: 2rem;
          font-weight: bold;
          color: #3182ce;
        }
        
        .stat-label {
          font-size: 14px;
          color: #4a5568;
          margin-top: 5px;
        }
        
        @media (max-width: 768px) {
          .experiments-grid {
            grid-template-columns: 1fr;
          }
          
          .experiment-header {
            flex-direction: column;
            align-items: flex-start;
            gap: 8px;
          }
          
          .page-title {
            font-size: 2rem;
          }
          
          .search-section {
            flex-direction: column;
            align-items: stretch;
          }
          
          .search-input {
            min-width: auto;
          }
        }
      `}</style>
    </div>
  );
}

// This function gets called at build time on server-side.
export async function getStaticProps() {
  return {
    props: {
      initialExperiments: mockExperiments
    },
  }
}

export default function Experiments({ initialExperiments }) {
  const [experiments, setExperiments] = useState(initialExperiments);
  const [loading, setLoading] = useState(false);

  return (
    <Layout title="Experiments">
      {isConvexEnabled ? (
        <ConvexExperimentsView />
      ) : (
        <StandardExperimentsView experiments={experiments} />
      )}
    </Layout>
  );
}