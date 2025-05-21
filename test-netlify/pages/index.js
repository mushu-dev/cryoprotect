import React from 'react';

export default function Home() {
  return (
    <div style={{ padding: '20px', maxWidth: '800px', margin: '0 auto' }}>
      <h1>CryoProtect on Netlify</h1>
      <p>This is a test deployment to verify Netlify can deploy our Next.js app.</p>
      
      <div style={{ marginTop: '20px', padding: '15px', backgroundColor: '#f0f0f0', borderRadius: '5px' }}>
        <h2>API Connection Test</h2>
        <p>Click the button below to test the connection to the Heroku backend.</p>
        <button 
          onClick={async () => {
            try {
              const response = await fetch('https://cryoprotect-8030e4025428.herokuapp.com/v1/health');
              if (response.ok) {
                alert('Successfully connected to the Heroku backend!');
              } else {
                alert(`Failed to connect to the backend. Status: ${response.status}`);
              }
            } catch (error) {
              alert(`Error connecting to backend: ${error.message}`);
            }
          }}
          style={{ 
            padding: '10px 15px', 
            backgroundColor: '#0070f3', 
            color: 'white', 
            border: 'none', 
            borderRadius: '4px',
            cursor: 'pointer'
          }}
        >
          Test API Connection
        </button>
      </div>
      
      <div style={{ marginTop: '20px' }}>
        <h2>Environment Variables</h2>
        <ul>
          <li>NEXT_PUBLIC_API_URL: {process.env.NEXT_PUBLIC_API_URL || 'Not set'}</li>
          <li>NEXT_PUBLIC_ENVIRONMENT: {process.env.NEXT_PUBLIC_ENVIRONMENT || 'Not set'}</li>
        </ul>
      </div>
    </div>
  );
}