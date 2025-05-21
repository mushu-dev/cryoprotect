import React from 'react';
import Head from 'next/head';

export default function Home() {
  return (
    <div style={{ 
      display: 'flex', 
      flexDirection: 'column', 
      alignItems: 'center', 
      justifyContent: 'center', 
      minHeight: '100vh',
      fontFamily: 'Arial, sans-serif',
      padding: '20px',
      textAlign: 'center'
    }}>
      <Head>
        <title>CryoProtect</title>
        <meta name="description" content="Cryoprotectant research and analysis platform" />
        <link rel="icon" href="/favicon.ico" />
      </Head>

      <h1 style={{ color: '#0070f3', marginBottom: '20px' }}>
        Welcome to CryoProtect
      </h1>

      <p style={{ fontSize: '1.2rem', maxWidth: '600px', lineHeight: 1.6 }}>
        CryoProtect is a platform for cryoprotectant research, analysis, and experimental tracking. 
        We're currently working on enhancements to our platform to better serve the scientific community.
      </p>

      <div style={{ 
        marginTop: '40px',
        padding: '20px',
        border: '1px solid #eaeaea',
        borderRadius: '10px',
        backgroundColor: '#f9f9f9',
        maxWidth: '600px'
      }}>
        <h2 style={{ color: '#0070f3', marginBottom: '10px' }}>Coming Soon</h2>
        <ul style={{ 
          textAlign: 'left', 
          paddingLeft: '20px',
          lineHeight: 1.6
        }}>
          <li>Enhanced molecular visualization</li>
          <li>Improved experimental data tracking</li>
          <li>Protocol and experiment sharing</li>
          <li>Advanced analysis tools</li>
          <li>Integration with laboratory equipment</li>
        </ul>
      </div>

      <footer style={{
        marginTop: '60px',
        color: '#666',
        fontSize: '0.9rem'
      }}>
        © {new Date().getFullYear()} CryoProtect - All rights reserved
      </footer>
    </div>
  );
}
