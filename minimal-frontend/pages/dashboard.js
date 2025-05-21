import React, { useEffect, useState } from 'react';
import { useRouter } from 'next/router';
import Link from 'next/link';
import Layout from '../components/Layout';
import SearchBar from '../components/SearchBar';
import { getMolecules, getMixtures, checkConvexConnection } from '../utils/api';

export default function Dashboard() {
  const router = useRouter();
  const [stats, setStats] = useState({
    molecules: { count: 0, loading: true, error: null },
    mixtures: { count: 0, loading: true, error: null },
    experiments: { count: 0, loading: true, error: null },
    protocols: { count: 0, loading: true, error: null },
    convexStatus: 'checking'
  });

  // Fetch dashboard data
  useEffect(() => {
    const fetchStats = async () => {
      try {
        // Check Convex connection status
        const isConvexConnected = await checkConvexConnection();
        setStats(prev => ({ ...prev, convexStatus: isConvexConnected ? 'online' : 'offline' }));
        
        // Get molecule count
        try {
          const molecules = await getMolecules();
          setStats(prev => ({ 
            ...prev, 
            molecules: { 
              count: molecules.length, 
              loading: false, 
              error: null 
            }
          }));
        } catch (error) {
          setStats(prev => ({ 
            ...prev, 
            molecules: { 
              count: 0, 
              loading: false, 
              error: error.message 
            }
          }));
        }
        
        // Get mixture count
        try {
          const mixtures = await getMixtures();
          setStats(prev => ({ 
            ...prev, 
            mixtures: { 
              count: mixtures.length, 
              loading: false, 
              error: null 
            }
          }));
        } catch (error) {
          setStats(prev => ({ 
            ...prev, 
            mixtures: { 
              count: 0, 
              loading: false, 
              error: error.message 
            }
          }));
        }
        
        // For demonstration purposes, we'll just set some demo values for experiments and protocols
        setStats(prev => ({
          ...prev,
          experiments: {
            count: 12,
            loading: false,
            error: null
          },
          protocols: {
            count: 5,
            loading: false,
            error: null
          }
        }));
      } catch (error) {
        console.error('Error fetching dashboard stats:', error);
      }
    };
    
    fetchStats();
  }, []);

  return (
    <Layout title="Dashboard">
      <div className="page-header">
        <h1 className="title">Dashboard</h1>
        <p className="description">
          Overview of your CryoProtect data and activities
        </p>
        
        <div className="dashboard-search">
          <SearchBar 
            placeholder="Search molecules, mixtures, and more..." 
            onSearch={(term, filter) => {
              router.push({
                pathname: '/search',
                query: { q: term, filter: filter || 'all' }
              });
            }}
          />
        </div>
      </div>
      
      {/* Stats Cards */}
      <div className="stats-grid">
        <div className="stat-card">
          <h2>Molecules</h2>
          <div className="stat-value">
            {stats.molecules.loading ? (
              <span className="loading">Loading...</span>
            ) : stats.molecules.error ? (
              <span className="error">Error: {stats.molecules.error}</span>
            ) : (
              <span>{stats.molecules.count}</span>
            )}
          </div>
          <Link href="/molecules" className="stat-link">
            View All Molecules
          </Link>
        </div>
        
        <div className="stat-card">
          <h2>Mixtures</h2>
          <div className="stat-value">
            {stats.mixtures.loading ? (
              <span className="loading">Loading...</span>
            ) : stats.mixtures.error ? (
              <span className="error">Error: {stats.mixtures.error}</span>
            ) : (
              <span>{stats.mixtures.count}</span>
            )}
          </div>
          <Link href="/mixtures" className="stat-link">
            View All Mixtures
          </Link>
        </div>
        
        <div className="stat-card">
          <h2>Experiments</h2>
          <div className="stat-value">
            {stats.experiments.loading ? (
              <span className="loading">Loading...</span>
            ) : stats.experiments.error ? (
              <span className="error">Error: {stats.experiments.error}</span>
            ) : (
              <span>{stats.experiments.count}</span>
            )}
          </div>
          <Link href="/experiments" className="stat-link">
            View All Experiments
          </Link>
        </div>
        
        <div className="stat-card">
          <h2>Protocols</h2>
          <div className="stat-value">
            {stats.protocols.loading ? (
              <span className="loading">Loading...</span>
            ) : stats.protocols.error ? (
              <span className="error">Error: {stats.protocols.error}</span>
            ) : (
              <span>{stats.protocols.count}</span>
            )}
          </div>
          <Link href="/protocols" className="stat-link">
            View All Protocols
          </Link>
        </div>
      </div>
      
      {/* Backend Status */}
      <div className="card">
        <h2 className="subtitle">Backend Status</h2>
        <div className="backend-status">
          <div className="status-item">
            <span className="status-label">Heroku API:</span>
            <span className={`status-value status-online`}>Online</span>
          </div>
          
          <div className="status-item">
            <span className="status-label">Fly.io RDKit:</span>
            <span className={`status-value status-online`}>Online</span>
          </div>
          
          <div className="status-item">
            <span className="status-label">Convex:</span>
            <span className={`status-value status-${stats.convexStatus}`}>
              {stats.convexStatus === 'checking' ? 'Checking...' : 
               stats.convexStatus === 'online' ? 'Online' : 'Offline'}
            </span>
          </div>
          
          <Link href="/connections" className="details-link">
            View Detailed Connection Status
          </Link>
        </div>
      </div>
      
      {/* Recent Activity */}
      <div className="card">
        <h2 className="subtitle">Recent Activity</h2>
        <div className="activity-list">
          <div className="activity-item">
            <div className="activity-icon">🧪</div>
            <div className="activity-content">
              <div className="activity-title">New Experiment Added</div>
              <div className="activity-description">Glycerol Freezing Point Experiment</div>
              <div className="activity-time">2 hours ago</div>
            </div>
          </div>
          
          <div className="activity-item">
            <div className="activity-icon">🧬</div>
            <div className="activity-content">
              <div className="activity-title">Mixture Updated</div>
              <div className="activity-description">VS55 Composition Modified</div>
              <div className="activity-time">Yesterday</div>
            </div>
          </div>
          
          <div className="activity-item">
            <div className="activity-icon">📋</div>
            <div className="activity-content">
              <div className="activity-title">Protocol Added</div>
              <div className="activity-description">Cell Line Freezing Protocol</div>
              <div className="activity-time">3 days ago</div>
            </div>
          </div>
        </div>
      </div>
      
      {/* Quick Actions */}
      <div className="card">
        <h2 className="subtitle">Quick Actions</h2>
        <div className="quick-actions">
          <Link href="/molecules" className="quick-action-button">
            Browse Molecules
          </Link>
          <Link href="/mixtures" className="quick-action-button">
            Browse Mixtures
          </Link>
          <Link href="/convex" className="quick-action-button">
            Convex Demo
          </Link>
          <Link href="/molecule-viewer" className="quick-action-button">
            Molecule Viewer
          </Link>
        </div>
      </div>
      
      <style jsx>{`
        .stats-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
          gap: 20px;
          margin-bottom: 30px;
        }
        
        .stat-card {
          padding: 20px;
          background-color: #f9f9f9;
          border-radius: 8px;
          box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
          text-align: center;
        }
        
        .stat-card h2 {
          margin-top: 0;
          margin-bottom: 15px;
          font-size: 18px;
          color: #333;
        }
        
        .stat-value {
          font-size: 36px;
          font-weight: bold;
          color: #3182ce;
          margin-bottom: 15px;
        }
        
        .stat-link, .details-link {
          display: inline-block;
          margin-top: 10px;
          font-size: 14px;
          color: #3182ce;
          text-decoration: none;
        }
        
        .stat-link:hover, .details-link:hover {
          text-decoration: underline;
        }
        
        .loading {
          font-size: 18px;
          color: #718096;
        }
        
        .error {
          font-size: 14px;
          color: #e53e3e;
        }
        
        .backend-status {
          margin-top: 15px;
        }
        
        .status-item {
          display: flex;
          justify-content: space-between;
          padding: 10px;
          border-bottom: 1px solid #eee;
        }
        
        .status-label {
          font-weight: 500;
        }
        
        .status-value {
          padding: 2px 8px;
          border-radius: 4px;
          font-size: 14px;
        }
        
        .status-online {
          background-color: #9ae6b4;
          color: #22543d;
        }
        
        .status-offline {
          background-color: #feb2b2;
          color: #742a2a;
        }
        
        .status-checking {
          background-color: #e9d8fd;
          color: #44337a;
        }
        
        .activity-list {
          margin-top: 15px;
        }
        
        .activity-item {
          display: flex;
          padding: 15px 10px;
          border-bottom: 1px solid #eee;
          transition: background-color 0.2s ease;
        }
        
        .activity-item:hover {
          background-color: #f9f9f9;
        }
        
        .activity-icon {
          font-size: 24px;
          margin-right: 15px;
        }
        
        .activity-content {
          flex: 1;
        }
        
        .activity-title {
          font-weight: 500;
          margin-bottom: 5px;
        }
        
        .activity-description {
          color: #4a5568;
          margin-bottom: 5px;
        }
        
        .activity-time {
          font-size: 12px;
          color: #718096;
        }
        
        .quick-actions {
          display: flex;
          flex-wrap: wrap;
          gap: 10px;
          margin-top: 15px;
          justify-content: center;
        }
        
        .quick-action-button {
          padding: 12px 18px;
          background-color: #3182ce;
          color: white;
          border-radius: 6px;
          text-decoration: none;
          font-weight: 500;
          font-size: 14px;
          transition: all 0.2s;
          box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
          text-align: center;
          min-width: 120px;
        }
        
        .quick-action-button:hover {
          background-color: #2c5282;
          transform: translateY(-2px);
          box-shadow: 0 4px 8px rgba(0, 0, 0, 0.15);
        }
        
        .dashboard-search {
          margin-top: 20px;
          max-width: 600px;
          width: 100%;
        }
      `}</style>
    </Layout>
  );
}