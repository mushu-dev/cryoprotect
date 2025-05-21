exports.handler = async function(event, context) {
  return {
    statusCode: 200,
    body: JSON.stringify({
      message: "CryoProtect API is running!",
      environment: process.env.NEXT_PUBLIC_ENVIRONMENT,
      timestamp: new Date().toISOString()
    })
  };
};