export default function Home() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '800px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2.5rem', marginBottom: '20px' }}>CryoProtect</h1>
      <p style={{ marginBottom: '30px' }}>Cryoprotectant Analysis Platform</p>
      <div style={{ display: 'flex', justifyContent: 'center', gap: '20px', flexWrap: 'wrap' }}>
        <a 
          href="/molecules"
          style={{ 
            display: 'inline-block', 
            padding: '10px 20px', 
            backgroundColor: '#0070f3', 
            color: 'white', 
            borderRadius: '5px', 
            textDecoration: 'none' 
          }}
        >
          Explore Molecules
        </a>
        <a 
          href="/mixtures"
          style={{ 
            display: 'inline-block', 
            padding: '10px 20px', 
            backgroundColor: '#0070f3', 
            color: 'white', 
            borderRadius: '5px', 
            textDecoration: 'none' 
          }}
        >
          View Mixtures
        </a>
      </div>
      <p style={{ marginTop: '40px', fontSize: '0.9rem', color: '#666' }}>
        This is a static fallback page. Please enable JavaScript for the full experience.
      </p>
    </div>
  )
}
