/**
 * ConvexExperimentsView Component
 * 
 * This component displays experiments from Convex with real-time updates.
 * It provides search, filtering, and detailed experiment cards.
 */
import React, { useState, useEffect } from 'react';
import { useExperiments } from '../convex/hooks';

const ConvexExperimentsView = () => {
  const [searchTerm, setSearchTerm] = useState('');
  const [statusFilter, setStatusFilter] = useState('all');
  const [typeFilter, setTypeFilter] = useState('all');
  
  // Query experiments from Convex with real-time updates
  const experiments = useExperiments(100); // Fetch up to 100 experiments
  
  // Filter experiments based on search and filters
  const filteredExperiments = React.useMemo(() => {
    if (!experiments) return [];
    
    return experiments.filter(experiment => {
      // Search filter
      const matchesSearch = !searchTerm || (
        experiment.name?.toLowerCase().includes(searchTerm.toLowerCase()) ||
        experiment.description?.toLowerCase().includes(searchTerm.toLowerCase()) ||
        experiment.researcher?.toLowerCase().includes(searchTerm.toLowerCase())
      );
      
      // Status filter
      const matchesStatus = statusFilter === 'all' || experiment.status === statusFilter;
      
      // Type filter
      const matchesType = typeFilter === 'all' || experiment.experimentType === typeFilter;
      
      return matchesSearch && matchesStatus && matchesType;
    });
  }, [experiments, searchTerm, statusFilter, typeFilter]);
  
  // Loading state
  if (experiments === undefined) {
    return (
      <div className="convex-experiments-view">
        <style jsx>{`
          .convex-experiments-view {
            padding: 2rem;
            text-align: center;
          }
          .loading {
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 1rem;
            color: #666;
          }
          .spinner {
            width: 40px;
            height: 40px;
            border: 3px solid #f3f3f3;
            border-top: 3px solid #3498db;
            border-radius: 50%;
            animation: spin 1s linear infinite;
          }
          @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
          }
        `}</style>
        <div className="loading">
          <div className="spinner"></div>
          <p>Loading experiments from Convex...</p>
        </div>
      </div>
    );
  }
  
  // Calculate statistics
  const stats = React.useMemo(() => {
    if (!experiments) return { total: 0, completed: 0, inProgress: 0, planned: 0 };
    
    return {
      total: experiments.length,
      completed: experiments.filter(exp => exp.status === 'completed').length,
      inProgress: experiments.filter(exp => exp.status === 'in_progress').length,
      planned: experiments.filter(exp => exp.status === 'planned').length
    };
  }, [experiments]);
  
  return (
    <div className="convex-experiments-view">
      <style jsx>{`
        .convex-experiments-view {
          padding: 2rem;
          max-width: 1200px;
          margin: 0 auto;
        }
        
        .header {
          margin-bottom: 2rem;
        }
        
        .header h1 {
          font-size: 2rem;
          font-weight: 700;
          color: #2c3e50;
          margin: 0 0 0.5rem 0;
        }
        
        .header .subtitle {
          color: #7f8c8d;
          font-size: 1.1rem;
          margin-bottom: 1.5rem;
        }
        
        .stats-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
          gap: 1rem;
          margin-bottom: 2rem;
        }
        
        .stat-card {
          background: white;
          padding: 1.5rem;
          border-radius: 12px;
          box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
          text-align: center;
          border-left: 4px solid #3498db;
        }
        
        .stat-card.completed { border-left-color: #27ae60; }
        .stat-card.in-progress { border-left-color: #f39c12; }
        .stat-card.planned { border-left-color: #9b59b6; }
        
        .stat-number {
          font-size: 2rem;
          font-weight: 700;
          color: #2c3e50;
          margin-bottom: 0.5rem;
        }
        
        .stat-label {
          color: #7f8c8d;
          font-size: 0.9rem;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }
        
        .controls {
          display: flex;
          flex-wrap: wrap;
          gap: 1rem;
          margin-bottom: 2rem;
          align-items: center;
        }
        
        .search-box {
          flex: 1;
          min-width: 250px;
          position: relative;
        }
        
        .search-box input {
          width: 100%;
          padding: 0.75rem 1rem;
          border: 2px solid #e1e8ed;
          border-radius: 25px;
          font-size: 1rem;
          outline: none;
          transition: border-color 0.2s ease;
        }
        
        .search-box input:focus {
          border-color: #3498db;
        }
        
        .filter-group {
          display: flex;
          gap: 1rem;
          align-items: center;
        }
        
        .filter-group select {
          padding: 0.5rem 1rem;
          border: 2px solid #e1e8ed;
          border-radius: 8px;
          font-size: 0.9rem;
          outline: none;
          background: white;
          cursor: pointer;
        }
        
        .filter-group select:focus {
          border-color: #3498db;
        }
        
        .experiments-grid {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
          gap: 1.5rem;
        }
        
        .experiment-card {
          background: white;
          border-radius: 12px;
          box-shadow: 0 4px 15px rgba(0, 0, 0, 0.1);
          padding: 1.5rem;
          transition: transform 0.2s ease, box-shadow 0.2s ease;
          cursor: pointer;
        }
        
        .experiment-card:hover {
          transform: translateY(-2px);
          box-shadow: 0 8px 25px rgba(0, 0, 0, 0.15);
        }
        
        .experiment-header {
          display: flex;
          justify-content: space-between;
          align-items: flex-start;
          margin-bottom: 1rem;
        }
        
        .experiment-title {
          font-size: 1.2rem;
          font-weight: 600;
          color: #2c3e50;
          margin: 0;
          flex: 1;
        }
        
        .status-badge {
          padding: 0.25rem 0.75rem;
          border-radius: 20px;
          font-size: 0.8rem;
          font-weight: 500;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          white-space: nowrap;
          margin-left: 1rem;
        }
        
        .status-planned { background: #e3f2fd; color: #1976d2; }
        .status-in_progress { background: #fff3e0; color: #f57c00; }
        .status-completed { background: #e8f5e8; color: #388e3c; }
        .status-aborted { background: #ffebee; color: #d32f2f; }
        .status-failed { background: #ffebee; color: #d32f2f; }
        
        .experiment-description {
          color: #7f8c8d;
          margin-bottom: 1rem;
          line-height: 1.4;
        }
        
        .experiment-details {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 0.5rem;
          margin-bottom: 1rem;
          font-size: 0.9rem;
        }
        
        .detail-item {
          display: flex;
          justify-content: space-between;
        }
        
        .detail-label {
          color: #7f8c8d;
        }
        
        .detail-value {
          color: #2c3e50;
          font-weight: 500;
        }
        
        .experiment-footer {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-top: 1rem;
          padding-top: 1rem;
          border-top: 1px solid #ecf0f1;
        }
        
        .researcher {
          color: #7f8c8d;
          font-size: 0.9rem;
        }
        
        .experiment-date {
          color: #7f8c8d;
          font-size: 0.9rem;
        }
        
        .no-experiments {
          text-align: center;
          padding: 3rem;
          color: #7f8c8d;
        }
        
        .no-experiments h3 {
          margin-bottom: 1rem;
          color: #95a5a6;
        }
        
        .convex-badge {
          display: inline-flex;
          align-items: center;
          gap: 0.5rem;
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          color: white;
          padding: 0.5rem 1rem;
          border-radius: 20px;
          font-size: 0.9rem;
          font-weight: 500;
          margin-bottom: 1rem;
        }
        
        .realtime-indicator {
          width: 8px;
          height: 8px;
          background: #27ae60;
          border-radius: 50%;
          animation: pulse 2s infinite;
        }
        
        @keyframes pulse {
          0% { opacity: 1; }
          50% { opacity: 0.5; }
          100% { opacity: 1; }
        }
      `}</style>
      
      <div className="header">
        <div className="convex-badge">
          <div className="realtime-indicator"></div>
          Real-time Convex Data
        </div>
        <h1>Experiments Database</h1>
        <p className="subtitle">Real-time experiment tracking and management</p>
      </div>
      
      <div className="stats-grid">
        <div className="stat-card">
          <div className="stat-number">{stats.total}</div>
          <div className="stat-label">Total Experiments</div>
        </div>
        <div className="stat-card completed">
          <div className="stat-number">{stats.completed}</div>
          <div className="stat-label">Completed</div>
        </div>
        <div className="stat-card in-progress">
          <div className="stat-number">{stats.inProgress}</div>
          <div className="stat-label">In Progress</div>
        </div>
        <div className="stat-card planned">
          <div className="stat-number">{stats.planned}</div>
          <div className="stat-label">Planned</div>
        </div>
      </div>
      
      <div className="controls">
        <div className="search-box">
          <input
            type="text"
            placeholder="Search experiments by name, description, or researcher..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
          />
        </div>
        
        <div className="filter-group">
          <select
            value={statusFilter}
            onChange={(e) => setStatusFilter(e.target.value)}
          >
            <option value="all">All Status</option>
            <option value="planned">Planned</option>
            <option value="in_progress">In Progress</option>
            <option value="completed">Completed</option>
            <option value="aborted">Aborted</option>
            <option value="failed">Failed</option>
          </select>
          
          <select
            value={typeFilter}
            onChange={(e) => setTypeFilter(e.target.value)}
          >
            <option value="all">All Types</option>
            <option value="cryopreservation">Cryopreservation</option>
            <option value="vitrification">Vitrification</option>
            <option value="freezing">Freezing</option>
            <option value="thawing">Thawing</option>
            <option value="storage">Storage</option>
            <option value="viability_study">Viability Study</option>
            <option value="recovery_study">Recovery Study</option>
            <option value="other">Other</option>
          </select>
        </div>
      </div>
      
      {filteredExperiments.length === 0 ? (
        <div className="no-experiments">
          <h3>No experiments found</h3>
          <p>
            {experiments?.length === 0 
              ? "No experiments in the database yet."
              : "No experiments match your search criteria."
            }
          </p>
        </div>
      ) : (
        <div className="experiments-grid">
          {filteredExperiments.map((experiment) => (
            <div 
              key={experiment._id} 
              className="experiment-card"
              onClick={() => {
                // Handle experiment card click - could navigate to detail page
                console.log('Clicked experiment:', experiment);
              }}
            >
              <div className="experiment-header">
                <h3 className="experiment-title">{experiment.name}</h3>
                <div className={`status-badge status-${experiment.status || 'planned'}`}>
                  {(experiment.status || 'planned').replace('_', ' ')}
                </div>
              </div>
              
              {experiment.description && (
                <p className="experiment-description">
                  {experiment.description.length > 150 
                    ? `${experiment.description.substring(0, 150)}...` 
                    : experiment.description
                  }
                </p>
              )}
              
              <div className="experiment-details">
                <div className="detail-item">
                  <span className="detail-label">Type:</span>
                  <span className="detail-value">
                    {experiment.experimentType || 'Not specified'}
                  </span>
                </div>
                <div className="detail-item">
                  <span className="detail-label">Tissue:</span>
                  <span className="detail-value">
                    {experiment.tissueType || 'Not specified'}
                  </span>
                </div>
                <div className="detail-item">
                  <span className="detail-label">Cryoprotectant:</span>
                  <span className="detail-value">
                    {experiment.cryoprotectant || 'Not specified'}
                  </span>
                </div>
                <div className="detail-item">
                  <span className="detail-label">Concentration:</span>
                  <span className="detail-value">
                    {experiment.concentration || 'Not specified'}
                  </span>
                </div>
              </div>
              
              <div className="experiment-footer">
                <div className="researcher">
                  {experiment.researcher && `Researcher: ${experiment.researcher}`}
                </div>
                <div className="experiment-date">
                  {experiment.startDate && 
                    `Started: ${new Date(experiment.startDate).toLocaleDateString()}`
                  }
                </div>
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
};

export default ConvexExperimentsView;