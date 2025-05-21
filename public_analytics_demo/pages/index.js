import { Analytics } from '@vercel/analytics/react'
import { SpeedInsights } from '@vercel/speed-insights/next'

export default function Home() {
  return (
    <div style={{ 
      display: 'flex', 
      flexDirection: 'column', 
      alignItems: 'center', 
      justifyContent: 'center', 
      height: '100vh',
      padding: '0 2rem',
      textAlign: 'center'
    }}>
      <h1 style={{ fontSize: '2.5rem', marginBottom: '1rem' }}>CryoProtect Analytics Demo</h1>
      <p style={{ fontSize: '1.2rem', marginBottom: '2rem', maxWidth: '600px' }}>
        This page demonstrates Vercel Analytics and Speed Insights integration.
        It's a minimal deployment to test these features.
      </p>
      <div style={{ 
        background: '#f4f4f4', 
        padding: '1.5rem', 
        borderRadius: '8px',
        maxWidth: '600px'
      }}>
        <h2 style={{ marginTop: 0 }}>Features Enabled:</h2>
        <ul style={{ textAlign: 'left' }}>
          <li>Vercel Analytics</li>
          <li>Vercel Speed Insights</li>
        </ul>
        <p>Visit your Vercel dashboard to see the analytics data.</p>
      </div>
      
      {/* Include Analytics and SpeedInsights */}
      <Analytics />
      <SpeedInsights />
    </div>
  )
}
