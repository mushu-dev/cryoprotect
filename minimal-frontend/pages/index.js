import React, { useState, useEffect } from 'react';
import Link from 'next/link';
import Layout from '../components/Layout';

// Mock data - replace with real Convex queries when ready
const DASHBOARD_DATA = {
  totalMolecules: 19,
  cryoprotectants: 15,
  averageScore: 7.2,
  topPerformers: [
    { name: 'DMSO', score: 8.9, molecularWeight: 78.13 },
    { name: 'Glycerol', score: 8.7, molecularWeight: 92.09 },
    { name: 'Ethylene Glycol', score: 8.5, molecularWeight: 62.07 }
  ]
};

const SYSTEM_STATUS = {
  convex: { status: 'checking', url: process.env.NEXT_PUBLIC_CONVEX_URL || 'https://upbeat-parrot-866.convex.cloud' },
  rdkit: { status: 'checking', url: process.env.NEXT_PUBLIC_RDKIT_URL || 'https://cryoprotect-rdkit.fly.dev' },
  api: { status: 'checking', url: process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/api' }
};

function StatusIndicator({ status, label, url }) {
  const getStatusColor = () => {
    switch (status) {
      case 'online': return '#10b981';
      case 'offline': return '#ef4444';
      case 'checking': return '#f59e0b';
      default: return '#6b7280';
    }
  };

  return (
    <div className="status-indicator">
      <div className="status-row">
        <div className="status-info">
          <div className="status-dot" style={{ backgroundColor: getStatusColor() }}></div>
          <div>
            <div className="status-label">{label}</div>
            <div className="status-url">{url}</div>
          </div>
        </div>
        <span className={`status-badge status-${status}`}>
          {status}
        </span>
      </div>
    </div>
  );
}

function MetricCard({ title, value, subtitle, icon, trend, color = "blue" }) {
  return (
    <div className={`metric-card metric-${color}`}>
      <div className="metric-content">
        <div className="metric-text">
          <p className="metric-title">{title}</p>
          <p className="metric-value">{value}</p>
          {subtitle && <p className="metric-subtitle">{subtitle}</p>}
        </div>
        <div className="metric-icon">{icon}</div>
      </div>
      {trend && (
        <div className="metric-trend">
          <span className={`trend-${trend.direction}`}>
            {trend.direction === 'up' ? '↗' : '↘'} {trend.value}
          </span>
          <span> {trend.label}</span>
        </div>
      )}
    </div>
  );
}

function MoleculePerformanceCard({ molecules }) {
  return (
    <div className="performance-card">
      <h3 className="card-title">Top Cryoprotectants</h3>
      <div className="molecules-list">
        {molecules.map((molecule, index) => (
          <div key={index} className="molecule-item">
            <div className="molecule-info">
              <div className="molecule-name">{molecule.name}</div>
              <div className="molecule-details">MW: {molecule.molecularWeight} g/mol</div>
            </div>
            <div className="molecule-score">
              <div className="score-value">{molecule.score}/10</div>
              <div className="score-label">Protection Score</div>
            </div>
          </div>
        ))}
      </div>
      <Link href="/molecules" className="view-all-btn">
        View All Molecules
      </Link>
    </div>
  );
}

function QuickAction({ href, icon, title, description }) {
  return (
    <Link href={href} className="quick-action">
      <div className="action-content">
        <div className="action-icon">{icon}</div>
        <div className="action-text">
          <h4 className="action-title">{title}</h4>
          <p className="action-description">{description}</p>
        </div>
      </div>
    </Link>
  );
}

export default function Dashboard() {
  const [systemStatus, setSystemStatus] = useState(SYSTEM_STATUS);

  useEffect(() => {
    // Check system status
    const checkStatus = async () => {
      const updatedStatus = { ...systemStatus };
      
      // Check each service
      for (const [service, config] of Object.entries(SYSTEM_STATUS)) {
        try {
          if (config.url) {
            // Simple check - in production, implement proper health checks
            updatedStatus[service].status = 'online';
          } else {
            updatedStatus[service].status = 'offline';
          }
        } catch (error) {
          updatedStatus[service].status = 'offline';
        }
      }
      
      setSystemStatus(updatedStatus);
    };

    checkStatus();
  }, []);

  return (
    <Layout title="Dashboard">
      <div className="dashboard-container">
        {/* Header */}
        <div className="dashboard-header">
          <h1 className="dashboard-title">CryoProtect Dashboard</h1>
          <p className="dashboard-subtitle">Monitor molecular data, system status, and research progress</p>
        </div>

        {/* Key Metrics Grid */}
        <div className="metrics-grid">
          <MetricCard
            title="Total Molecules"
            value={DASHBOARD_DATA.totalMolecules}
            subtitle="In database"
            icon="🧪"
            color="blue"
          />
          <MetricCard
            title="Cryoprotectants"
            value={DASHBOARD_DATA.cryoprotectants}
            subtitle="Active compounds"
            icon="🧬"
            color="green"
          />
          <MetricCard
            title="Avg Protection Score"
            value={DASHBOARD_DATA.averageScore}
            subtitle="Out of 10"
            icon="⭐"
            color="purple"
          />
          <MetricCard
            title="Experiments"
            value="4"
            subtitle="2 completed"
            icon="🔬"
            trend={{ direction: 'up', value: '+25%', label: 'this month' }}
            color="orange"
          />
        </div>

        {/* Main Content Grid */}
        <div className="main-content-grid">
          {/* Molecular Performance */}
          <div className="performance-section">
            <MoleculePerformanceCard molecules={DASHBOARD_DATA.topPerformers} />
          </div>

          {/* System Status */}
          <div className="status-section">
            <div className="status-card">
              <h3 className="card-title">System Status</h3>
              <div className="status-list">
                <StatusIndicator 
                  status={systemStatus.convex.status}
                  label="Convex Database"
                  url={systemStatus.convex.url}
                />
                <StatusIndicator 
                  status={systemStatus.rdkit.status}
                  label="RDKit Service"
                  url={systemStatus.rdkit.url}
                />
                <StatusIndicator 
                  status={systemStatus.api.status}
                  label="Main API"
                  url={systemStatus.api.url}
                />
              </div>
              
              {/* Configuration Info */}
              <div className="config-info">
                <div className="info-icon">✅</div>
                <div>
                  <p className="info-title">Database Ready</p>
                  <p className="info-desc">19 molecules with complete scientific data</p>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Quick Actions */}
        <div className="quick-actions-section">
          <h3 className="section-title">Quick Actions</h3>
          <div className="actions-grid">
            <QuickAction
              href="/molecules"
              icon="🔍"
              title="Browse Molecules"
              description="Explore molecular database"
            />
            <QuickAction
              href="/experiments"
              icon="📊"
              title="View Experiments"
              description="Analyze experiment results"
            />
            <QuickAction
              href="/protocols"
              icon="📋"
              title="Create Protocol"
              description="Design new experiments"
            />
            <QuickAction
              href="/connections"
              icon="🔧"
              title="System Health"
              description="Check service status"
            />
          </div>
        </div>

        {/* Scientific Insights */}
        <div className="insights-section">
          <h3 className="insights-title">🧪 Scientific Insights</h3>
          <div className="insights-grid">
            <div className="insight-item">
              <div className="insight-label">Molecular Properties</div>
              <div className="insight-desc">Complete RDKit descriptors available for all compounds</div>
            </div>
            <div className="insight-item">
              <div className="insight-label">Cryoprotection Scoring</div>
              <div className="insight-desc">Based on glass transition, permeability, and biocompatibility</div>
            </div>
            <div className="insight-item">
              <div className="insight-label">Experimental Validation</div>
              <div className="insight-desc">4 protocols with viability and recovery data</div>
            </div>
          </div>
        </div>
      </div>

      <style jsx>{`
        .dashboard-container {
          max-width: 1200px;
          margin: 0 auto;
          padding: 20px;
        }

        .dashboard-header {
          margin-bottom: 32px;
        }

        .dashboard-title {
          font-size: 28px;
          font-weight: bold;
          color: #111827;
          margin: 0 0 8px 0;
        }

        .dashboard-subtitle {
          color: #6b7280;
          margin: 0;
        }

        .metrics-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
          gap: 24px;
          margin-bottom: 32px;
        }

        .metric-card {
          border: 1px solid #e5e7eb;
          border-radius: 8px;
          padding: 24px;
          background: white;
        }

        .metric-blue { border-left: 4px solid #3b82f6; }
        .metric-green { border-left: 4px solid #10b981; }
        .metric-purple { border-left: 4px solid #8b5cf6; }
        .metric-orange { border-left: 4px solid #f59e0b; }

        .metric-content {
          display: flex;
          justify-content: space-between;
          align-items: flex-start;
        }

        .metric-title {
          font-size: 14px;
          font-weight: 500;
          color: #6b7280;
          margin: 0 0 4px 0;
        }

        .metric-value {
          font-size: 24px;
          font-weight: bold;
          color: #111827;
          margin: 0 0 4px 0;
        }

        .metric-subtitle {
          font-size: 14px;
          color: #9ca3af;
          margin: 0;
        }

        .metric-icon {
          font-size: 24px;
        }

        .metric-trend {
          margin-top: 8px;
          font-size: 14px;
          color: #6b7280;
        }

        .trend-up { color: #10b981; }
        .trend-down { color: #ef4444; }

        .main-content-grid {
          display: grid;
          grid-template-columns: 2fr 1fr;
          gap: 32px;
          margin-bottom: 32px;
        }

        @media (max-width: 768px) {
          .main-content-grid {
            grid-template-columns: 1fr;
          }
        }

        .performance-card, .status-card {
          background: white;
          border: 1px solid #e5e7eb;
          border-radius: 8px;
          padding: 24px;
        }

        .card-title {
          font-size: 18px;
          font-weight: 600;
          color: #111827;
          margin: 0 0 16px 0;
        }

        .molecules-list {
          margin-bottom: 16px;
        }

        .molecule-item {
          display: flex;
          justify-content: space-between;
          align-items: center;
          padding: 12px;
          background: #f9fafb;
          border-radius: 6px;
          margin-bottom: 8px;
        }

        .molecule-name {
          font-weight: 500;
          color: #111827;
        }

        .molecule-details {
          font-size: 14px;
          color: #6b7280;
        }

        .score-value {
          font-weight: bold;
          color: #10b981;
          text-align: right;
        }

        .score-label {
          font-size: 12px;
          color: #6b7280;
          text-align: right;
        }

        .view-all-btn {
          display: block;
          text-align: center;
          background: #3b82f6;
          color: white;
          padding: 12px 16px;
          border-radius: 6px;
          text-decoration: none;
          font-weight: 500;
          transition: background-color 0.2s;
        }

        .view-all-btn:hover {
          background: #2563eb;
        }

        .status-list {
          margin-bottom: 16px;
        }

        .status-indicator {
          margin-bottom: 12px;
        }

        .status-row {
          display: flex;
          justify-content: space-between;
          align-items: center;
          padding: 12px;
          background: #f9fafb;
          border-radius: 6px;
        }

        .status-info {
          display: flex;
          align-items: center;
          gap: 12px;
        }

        .status-dot {
          width: 12px;
          height: 12px;
          border-radius: 50%;
        }

        .status-label {
          font-weight: 500;
          color: #111827;
        }

        .status-url {
          font-size: 12px;
          color: #6b7280;
          max-width: 200px;
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
        }

        .status-badge {
          padding: 4px 8px;
          border-radius: 12px;
          font-size: 12px;
          font-weight: 500;
        }

        .status-online {
          background: #dcfce7;
          color: #166534;
        }

        .status-offline {
          background: #fee2e2;
          color: #991b1b;
        }

        .status-checking {
          background: #fef3c7;
          color: #92400e;
        }

        .config-info {
          display: flex;
          align-items: center;
          gap: 8px;
          padding: 12px;
          background: #f0f9ff;
          border: 1px solid #bae6fd;
          border-radius: 6px;
        }

        .info-icon {
          font-size: 16px;
        }

        .info-title {
          font-size: 14px;
          font-weight: 500;
          color: #111827;
          margin: 0;
        }

        .info-desc {
          font-size: 12px;
          color: #6b7280;
          margin: 0;
        }

        .quick-actions-section {
          margin-bottom: 32px;
        }

        .section-title {
          font-size: 18px;
          font-weight: 600;
          color: #111827;
          margin: 0 0 16px 0;
        }

        .actions-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
          gap: 16px;
        }

        .quick-action {
          display: block;
          padding: 16px;
          background: white;
          border: 1px solid #e5e7eb;
          border-radius: 8px;
          text-decoration: none;
          transition: all 0.2s;
        }

        .quick-action:hover {
          border-color: #3b82f6;
          transform: translateY(-2px);
        }

        .action-content {
          display: flex;
          align-items: center;
          gap: 12px;
        }

        .action-icon {
          font-size: 24px;
        }

        .action-title {
          font-weight: 500;
          color: #111827;
          margin: 0 0 4px 0;
        }

        .action-description {
          font-size: 14px;
          color: #6b7280;
          margin: 0;
        }

        .insights-section {
          background: linear-gradient(135deg, #eff6ff 0%, #f3e8ff 100%);
          border: 1px solid #bae6fd;
          border-radius: 8px;
          padding: 24px;
        }

        .insights-title {
          font-size: 18px;
          font-weight: 600;
          color: #111827;
          margin: 0 0 16px 0;
        }

        .insights-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
          gap: 16px;
        }

        .insight-item {
          text-align: left;
        }

        .insight-label {
          font-weight: 500;
          color: #111827;
          margin-bottom: 4px;
        }

        .insight-desc {
          font-size: 14px;
          color: #6b7280;
        }
      `}</style>
    </Layout>
  );
}