import React from 'react';
import { Metadata } from 'next';

export const metadata: Metadata = {
  title: 'Experiments | CryoProtect',
  description: 'Design, track, and analyze cryopreservation experiments with detailed protocols',
};

export default function ExperimentsLayout({
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