export default function MixturePage() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>Mixture Details</h1>
      <p style={{ marginBottom: '30px' }}>Loading mixture information...</p>
      <p style={{ fontSize: '0.9rem', color: '#666' }}>
        Please enable JavaScript to view mixture details.
      </p>
      <a 
        href="/mixtures"
        style={{ 
          display: 'inline-block', 
          marginTop: '20px',
          padding: '10px 20px', 
          backgroundColor: '#0070f3', 
          color: 'white', 
          borderRadius: '5px', 
          textDecoration: 'none' 
        }}
      >
        Back to All Mixtures
      </a>
    </div>
  )
}

export function getStaticPaths() {
  return {
    paths: [
      { params: { id: 'placeholder' } },
      { params: { id: '1' } },
      { params: { id: '2' } },
      { params: { id: '3' } }
    ],
    fallback: true
  }
}

export function getStaticProps() {
  return {
    props: {}
  }
}
