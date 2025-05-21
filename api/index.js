// JavaScript fallback for the index route
// This is a backup in case the Python serverless function doesn't work

// Vercel serverless function
export default function handler(req, res) {
  // Set CORS headers
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Access-Control-Allow-Methods', 'GET, POST, OPTIONS');
  res.setHeader('Access-Control-Allow-Headers', 'Content-Type');
  
  // Handle OPTIONS method for CORS preflight
  if (req.method === 'OPTIONS') {
    res.status(200).end();
    return;
  }
  
  // Get current path
  const path = req.url || '/';
  
  // Health check endpoint
  if (path.includes('/api/v1/health')) {
    return res.status(200).json({
      status: 'ok',
      version: '1.0.0',
      timestamp: new Date().toISOString(),
      environment: 'vercel-production',
      backend: 'JavaScript fallback'
    });
  }
  
  // Default API response for unknown endpoints
  return res.status(404).json({
    error: 'Unknown API route',
    path: path,
    timestamp: new Date().toISOString(),
    backend: 'JavaScript fallback'
  });
}