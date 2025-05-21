// API proxy function to handle all /api/* routes
exports.handler = async function(event, context) {
  const path = event.path;
  console.log("API Request received at path:", path);
  
  // Handle specific API routes
  if (path.endsWith('/api/hello')) {
    return {
      statusCode: 200,
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        name: 'CryoProtect API',
        message: 'API is working!',
        timestamp: new Date().toISOString(),
        environment: process.env.NEXT_PUBLIC_ENVIRONMENT || 'production'
      })
    };
  }
  
  if (path.endsWith('/api/health') || path.endsWith('/health')) {
    return {
      statusCode: 200,
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        status: 'OK',
        timestamp: new Date().toISOString(),
        environment: process.env.NEXT_PUBLIC_ENVIRONMENT || 'production'
      })
    };
  }
  
  // Default response for unhandled API routes
  return {
    statusCode: 404,
    body: JSON.stringify({
      error: 'API route not found',
      path: path
    })
  };
};
