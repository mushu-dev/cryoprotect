export default function Custom404() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>404 - Page Not Found</h1>
      <p style={{ marginBottom: '30px' }}>The page you are looking for does not exist.</p>
      <a 
        href="/"
        style={{ 
          display: 'inline-block', 
          padding: '10px 20px', 
          backgroundColor: '#0070f3', 
          color: 'white', 
          borderRadius: '5px', 
          textDecoration: 'none' 
        }}
      >
        Back to Home
      </a>
    </div>
  )
}
