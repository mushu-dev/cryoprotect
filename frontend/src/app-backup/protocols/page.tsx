import React from 'react';
import type { Metadata } from 'next';
import ProtocolsList from '../../features/protocols/components/protocols-list';

export const metadata: Metadata = {
  title: 'Protocols | CryoProtect',
  description: 'Browse and manage standardized protocols for cryopreservation experiments',
};

export default function ProtocolsPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <ProtocolsList />
    </div>
  );
}