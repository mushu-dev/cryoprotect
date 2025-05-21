import { useState } from 'react';
import Head from 'next/head';

// This page is completely static with no imports from Convex to avoid build errors
export default function ConvexDemoPage() {
  const [showMessage, setShowMessage] = useState(false);

  return (
    <div>
      <Head>
        <title>Convex Integration | CryoProtect</title>
        <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        <meta name="description" content="Convex integration demo for CryoProtect" />
      </Head>

      <main style={{ 
        maxWidth: '1200px', 
        margin: '0 auto', 
        padding: '20px',
        fontFamily: 'sans-serif'
      }}>
        <div style={{ marginBottom: '20px' }}>
          <h1 style={{ fontSize: '2rem', marginBottom: '10px' }}>
            Convex Integration Demo
          </h1>
          <p style={{ fontSize: '1rem', color: '#666' }}>
            This is a placeholder page. Convex integration is disabled in this minimal deployment.
          </p>
        </div>
        
        <div style={{
          padding: '20px',
          border: '1px solid #eee',
          borderRadius: '8px',
          marginTop: '20px',
          backgroundColor: '#f9f9f9'
        }}>
          <h2 style={{ fontSize: '1.5rem', marginTop: 0 }}>Convex Features (Disabled)</h2>
          <p>Convex integration has been disabled for this deployment to simplify the build process.</p>
          <p>Please refer to the full application for Convex functionality.</p>
          
          <button 
            onClick={() => setShowMessage(!showMessage)}
            style={{
              padding: '8px 16px',
              backgroundColor: '#3182ce',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              marginTop: '10px',
              fontWeight: '500'
            }}
          >
            {showMessage ? 'Hide Message' : 'Show Message'}
          </button>
          
          {showMessage && (
            <div style={{ 
              marginTop: '15px', 
              padding: '10px', 
              backgroundColor: '#ebf8ff', 
              borderRadius: '4px',
              borderLeft: '4px solid #3182ce'
            }}>
              <p>Convex is a real-time database that allows for automatic synchronization between clients and the backend. 
              In the full application, this page displays synchronized data from our Convex database.</p>
            </div>
          )}
        </div>
        
        <div style={{ marginTop: '30px' }}>
          <a 
            href="/"
            style={{
              textDecoration: 'none',
              color: '#3182ce',
              fontWeight: '500'
            }}
          >
            ← Back to Home
          </a>
        </div>
      </main>
    </div>
  );
}