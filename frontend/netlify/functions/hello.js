exports.handler = async function(event, context) {
  return {
    statusCode: 200,
    body: JSON.stringify({
      message: "CryoProtect API is running!",
      timestamp: new Date().toISOString()
    })
  };
};