import ClientMixturePage from './client-page'

// Generate static params for static site generation
export function generateStaticParams() {
  return [
    { id: 'placeholder' },
    { id: '1' },  // Example mixture 1
    { id: '2' },  // Example mixture 2
    { id: '3' }   // Example mixture 3
  ]
}

// Mark this as a static page
export const dynamic = 'force-static';

export default function MixtureDetailPage({ params }: { params: { id: string } }) {
  // Server component that passes params to the client component
  return <ClientMixturePage id={params.id} />;
}