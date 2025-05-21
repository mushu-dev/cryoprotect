import { useState, useEffect } from 'react';

// Simplified ConnectionStatus component that doesn't try to import Convex
export default function ConnectionStatus() {
  const [herokuStatus, setHerokuStatus] = useState('checking');
  const [flyioStatus, setFlyioStatus] = useState('checking');

  // Check Heroku API status
  useEffect(() => {
    const checkHerokuApi = async () => {
      try {
        const response = await fetch('https://cryoprotect-8030e4025428.herokuapp.com/api/molecules', {
          method: 'GET',
          headers: {
            'Accept': 'application/json',
          },
          mode: 'cors',
          credentials: 'omit',
          timeout: 5000
        });
        
        if (response.ok) {
          setHerokuStatus('online');
        } else {
          setHerokuStatus('offline');
        }
      } catch (error) {
        console.error('Failed to connect to Heroku:', error);
        setHerokuStatus('offline');
      }
    };
    
    checkHerokuApi();
  }, []);

  // Check Fly.io RDKit service status
  useEffect(() => {
    const checkFlyIoService = async () => {
      try {
        const response = await fetch('https://cryoprotect-rdkit.fly.dev/health', {
          method: 'GET',
          headers: {
            'Accept': 'application/json',
          },
          mode: 'cors',
          credentials: 'omit',
          timeout: 5000
        });
        
        if (response.ok) {
          setFlyioStatus('online');
        } else {
          setFlyioStatus('offline');
        }
      } catch (error) {
        console.error('Failed to connect to Fly.io:', error);
        setFlyioStatus('offline');
      }
    };
    
    checkFlyIoService();
  }, []);

  return (
    <div className="connection-status">
      <h2>Backend Connection Status</h2>
      <div className="status-grid">
        <div className="status-item">
          <span className="status-label">Heroku API:</span>
          <span className={`status-value status-${herokuStatus}`}>
            {herokuStatus === 'checking' ? 'Checking...' : herokuStatus === 'online' ? 'Online' : 'Offline'}
          </span>
        </div>
        
        <div className="status-item">
          <span className="status-label">Fly.io RDKit:</span>
          <span className={`status-value status-${flyioStatus}`}>
            {flyioStatus === 'checking' ? 'Checking...' : flyioStatus === 'online' ? 'Online' : 'Offline'}
          </span>
        </div>
        
        <div className="status-item">
          <span className="status-label">Convex:</span>
          <span className="status-value status-not-configured">
            Not Configured
          </span>
        </div>
      </div>
    </div>
  );
}