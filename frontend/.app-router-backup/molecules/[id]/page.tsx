import ClientMoleculePage from './client-page'

// Generate static params for static site generation
export function generateStaticParams() {
  return [
    { id: 'placeholder' },
    { id: '962' },   // Glycerol
    { id: '176' },   // DMSO
    { id: '6276' },  // Ethylene Glycol 
    { id: '8857' }   // Propylene Glycol
  ]
}

// Mark this as a static page
export const dynamic = 'force-static';

export default function MoleculeDetailPage({ params }: { params: { id: string } }) {
  // Server component that passes params to the client component
  return <ClientMoleculePage id={params.id} />;
}