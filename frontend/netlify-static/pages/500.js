export default function Custom500() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>500 - Server Error</h1>
      <p style={{ marginBottom: '30px' }}>Something went wrong on our end. Please try again later.</p>
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
