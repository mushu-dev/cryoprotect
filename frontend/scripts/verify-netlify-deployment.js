// Script to verify the Netlify deployment connectivity to Heroku backend
const https = require('https');

const netlifyUrl = process.env.NETLIFY_URL || 'https://cryoprotect.netlify.app';
const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1';

console.log(`Verifying connectivity from ${netlifyUrl} to ${apiUrl}...`);

// Function to make an HTTP request and return a promise
function makeRequest(url) {
  return new Promise((resolve, reject) => {
    https.get(url, (res) => {
      let data = '';
      
      res.on('data', (chunk) => {
        data += chunk;
      });
      
      res.on('end', () => {
        resolve({
          statusCode: res.statusCode,
          headers: res.headers,
          data: data
        });
      });
      
    }).on('error', (err) => {
      reject(err);
    });
  });
}

// Check API health endpoint
async function checkApiHealth() {
  try {
    const healthEndpoint = `${apiUrl}/health`;
    console.log(`Checking API health at: ${healthEndpoint}`);
    
    const response = await makeRequest(healthEndpoint);
    
    console.log(`API Health Status Code: ${response.statusCode}`);
    console.log(`API Health Response: ${response.data.substring(0, 200)}...`);
    
    if (response.statusCode >= 200 && response.statusCode < 300) {
      console.log('✅ API health check successful');
      return true;
    } else {
      console.error('❌ API health check failed');
      return false;
    }
  } catch (error) {
    console.error('❌ API health check error:', error.message);
    return false;
  }
}

// Check Netlify deployment
async function checkNetlifyDeployment() {
  try {
    console.log(`Checking Netlify deployment at: ${netlifyUrl}`);
    
    const response = await makeRequest(netlifyUrl);
    
    console.log(`Netlify Status Code: ${response.statusCode}`);
    
    if (response.statusCode >= 200 && response.statusCode < 300) {
      console.log('✅ Netlify deployment check successful');
      return true;
    } else {
      console.error('❌ Netlify deployment check failed');
      return false;
    }
  } catch (error) {
    console.error('❌ Netlify deployment error:', error.message);
    return false;
  }
}

// Run the checks
async function runVerification() {
  const netlifyOk = await checkNetlifyDeployment();
  const apiOk = await checkApiHealth();
  
  if (netlifyOk && apiOk) {
    console.log('\n✅ Deployment verification successful!');
    process.exit(0);
  } else {
    console.error('\n❌ Deployment verification failed!');
    process.exit(1);
  }
}

runVerification();