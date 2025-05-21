export default function MoleculePage() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>Molecule Details</h1>
      <p style={{ marginBottom: '30px' }}>Loading molecule information...</p>
      <p style={{ fontSize: '0.9rem', color: '#666' }}>
        Please enable JavaScript to view molecule details.
      </p>
      <a 
        href="/molecules"
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
        Back to All Molecules
      </a>
    </div>
  )
}

export function getStaticPaths() {
  return {
    paths: [
      { params: { id: 'placeholder' } },
      { params: { id: '962' } },
      { params: { id: '176' } },
      { params: { id: '6276' } },
      { params: { id: '8857' } }
    ],
    fallback: true
  }
}

export function getStaticProps() {
  return {
    props: {}
  }
}
