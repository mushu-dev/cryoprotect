/**
 * Dynamic API route handler
 * 
 * This handler serves as a proxy for development and testing purposes
 * It supports simulating delays and errors for resilience testing
 */

export default async function handler(req, res) {
  // Parse parameters
  const { endpoint } = req.query;
  const delay = parseInt(req.query.delay) || 0;
  const shouldError = req.query.error === 'true';
  
  // Simulate processing delay if specified
  if (delay > 0) {
    await new Promise(resolve => setTimeout(resolve, delay));
  }
  
  // Simulate error if specified
  if (shouldError) {
    return res.status(500).json({ 
      error: 'Simulated API error',
      endpoint,
      timestamp: new Date().toISOString()
    });
  }
  
  // Mock endpoints for testing
  switch(endpoint) {
    case 'health':
      return res.status(200).json({ status: 'ok', timestamp: new Date().toISOString() });
    
    case 'test-cors':
      // Set CORS headers
      res.setHeader('Access-Control-Allow-Origin', '*');
      res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS,POST');
      res.setHeader('Access-Control-Allow-Headers', 'Content-Type');
      
      if (req.method === 'OPTIONS') {
        return res.status(200).end();
      }
      
      return res.status(200).json({ 
        cors: 'enabled',
        method: req.method,
        timestamp: new Date().toISOString()
      });
    
    case 'mock-molecules':
      return res.status(200).json({
        data: [
          { id: 1, name: 'DMSO', formula: 'C2H6OS', cid: '679' },
          { id: 2, name: 'Glycerol', formula: 'C3H8O3', cid: '753' },
          { id: 3, name: 'Ethylene glycol', formula: 'C2H6O2', cid: '174' }
        ],
        count: 3,
        timestamp: new Date().toISOString()
      });
    
    case 'mock-mixtures':
      return res.status(200).json({
        data: [
          { 
            id: 1, 
            name: 'Standard Vitrification Solution', 
            components: [
              { molecule_id: 1, name: 'DMSO', concentration: 10 },
              { molecule_id: 2, name: 'Glycerol', concentration: 10 }
            ]
          },
          { 
            id: 2, 
            name: 'Slow Freezing Solution', 
            components: [
              { molecule_id: 1, name: 'DMSO', concentration: 5 }
            ]
          }
        ],
        count: 2,
        timestamp: new Date().toISOString()
      });
      
    default:
      return res.status(404).json({ 
        error: 'Endpoint not found',
        endpoint,
        timestamp: new Date().toISOString()
      });
  }
}