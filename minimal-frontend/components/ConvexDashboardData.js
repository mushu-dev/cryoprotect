import React, { useState, useEffect } from 'react';

// Convex client setup
const CONVEX_URL = process.env.NEXT_PUBLIC_CONVEX_URL;

export function useConvexDashboardData() {
  const [data, setData] = useState({
    molecules: [],
    molecularProperties: [],
    cryoprotectantScores: [],
    loading: true,
    error: null
  });

  useEffect(() => {
    const fetchDashboardData = async () => {
      try {
        setData(prev => ({ ...prev, loading: true }));

        // Fetch molecules data from Convex
        const moleculesResponse = await fetch(`${CONVEX_URL}/api/query`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            path: 'molecules:list',
            args: {}
          })
        });

        if (!moleculesResponse.ok) {
          throw new Error(`Molecules fetch failed: ${moleculesResponse.status}`);
        }

        const molecules = await moleculesResponse.json();

        // Fetch molecular properties
        const propertiesResponse = await fetch(`${CONVEX_URL}/api/query`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            path: 'molecules:listMoleculesWithProperties',
            args: {}
          })
        });

        const properties = propertiesResponse.ok ? await propertiesResponse.json() : [];

        // Fetch cryoprotectant scores
        const scoresResponse = await fetch(`${CONVEX_URL}/api/query`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            path: 'molecules:listMoleculesWithScores',
            args: {}
          })
        });

        const scores = scoresResponse.ok ? await scoresResponse.json() : [];

        setData({
          molecules: molecules || [],
          molecularProperties: properties || [],
          cryoprotectantScores: scores || [],
          loading: false,
          error: null
        });

      } catch (error) {
        console.error('Error fetching dashboard data:', error);
        setData(prev => ({
          ...prev,
          loading: false,
          error: error.message
        }));
      }
    };

    if (CONVEX_URL) {
      fetchDashboardData();
    } else {
      setData(prev => ({
        ...prev,
        loading: false,
        error: 'Convex URL not configured'
      }));
    }
  }, []);

  return data;
}

export function processDashboardMetrics(data) {
  if (data.loading || data.error || !data.molecules) {
    return {
      totalMolecules: 0,
      cryoprotectants: 0,
      averageScore: 0,
      topPerformers: [],
      molecularClasses: {},
      loading: data.loading,
      error: data.error
    };
  }

  const molecules = data.molecules;
  const scores = data.cryoprotectantScores;
  const properties = data.molecularProperties;

  // Calculate metrics
  const totalMolecules = molecules.length;
  
  // Count cryoprotectants (molecules with scores)
  const cryoprotectants = scores.length;
  
  // Calculate average protection score
  const validScores = scores.filter(s => s.totalScore && s.totalScore > 0);
  const averageScore = validScores.length > 0 
    ? (validScores.reduce((sum, s) => sum + s.totalScore, 0) / validScores.length).toFixed(1)
    : 0;

  // Get top performers
  const topPerformers = scores
    .filter(s => s.totalScore && s.totalScore > 0)
    .sort((a, b) => b.totalScore - a.totalScore)
    .slice(0, 3)
    .map(score => {
      const molecule = molecules.find(m => m._id === score.moleculeId);
      const molProps = properties.find(p => p.moleculeId === score.moleculeId);
      
      return {
        name: molecule?.name || 'Unknown',
        score: parseFloat(score.totalScore).toFixed(1),
        molecularWeight: molProps?.molecularWeight || 0,
        smiles: molecule?.smiles || '',
        id: molecule?._id
      };
    });

  // Classify molecules by type
  const molecularClasses = molecules.reduce((acc, mol) => {
    const classification = mol.classification || 'Unclassified';
    acc[classification] = (acc[classification] || 0) + 1;
    return acc;
  }, {});

  return {
    totalMolecules,
    cryoprotectants,
    averageScore: parseFloat(averageScore),
    topPerformers,
    molecularClasses,
    loading: false,
    error: null
  };
}

export function ConvexConnectionStatus() {
  const [status, setStatus] = useState('checking');
  const [responseTime, setResponseTime] = useState(null);

  useEffect(() => {
    const checkConnection = async () => {
      if (!CONVEX_URL) {
        setStatus('offline');
        return;
      }

      try {
        const startTime = Date.now();
        
        const response = await fetch(`${CONVEX_URL}/api/query`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            path: 'molecules:list',
            args: {}
          }),
          timeout: 10000
        });

        const endTime = Date.now();
        const responseTimeMs = endTime - startTime;

        if (response.ok) {
          setStatus('online');
          setResponseTime(responseTimeMs);
        } else {
          setStatus('offline');
        }
      } catch (error) {
        console.error('Convex connection check failed:', error);
        setStatus('offline');
      }
    };

    checkConnection();
    
    // Check every 30 seconds
    const interval = setInterval(checkConnection, 30000);
    return () => clearInterval(interval);
  }, []);

  return {
    status,
    responseTime,
    url: CONVEX_URL
  };
}