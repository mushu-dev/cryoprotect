import { Suspense } from 'react';
import SettingsClientPage from './page.client';

export const dynamic = 'force-dynamic';

export default function SettingsPage() {
  return (
    <Suspense fallback={<div>Loading settings...</div>}>
      <SettingsClientPage />
    </Suspense>
  );
}