import { Suspense } from 'react';
import ProfileClientPage from './page.client';

export const dynamic = 'force-dynamic';

export default function ProfilePage() {
  return (
    <Suspense fallback={<div>Loading profile...</div>}>
      <ProfileClientPage />
    </Suspense>
  );
}