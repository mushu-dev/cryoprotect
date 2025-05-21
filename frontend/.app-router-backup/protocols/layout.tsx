import React from 'react';
import { Metadata } from 'next';

export const metadata: Metadata = {
  title: 'Protocols | CryoProtect',
  description: 'Browse and manage standardized protocols for cryopreservation experiments',
};

export default function ProtocolsLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <div>
      {children}
    </div>
  );
}