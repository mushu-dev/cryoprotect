import React from 'react';
import Link from 'next/link';
import Layout from '../components/Layout';

export default function NotFound() {
  return (
    <Layout title="Page Not Found">
      <div className="error-container">
        <h1 className="error-title">404</h1>
        <h2 className="error-subtitle">Page Not Found</h2>
        <p className="error-description">
          Sorry, the page you are looking for does not exist or has been moved.
        </p>
        <div className="error-actions">
          <Link href="/">
            <span className="button">Return to Home</span>
          </Link>
        </div>
      </div>
    </Layout>
  );
}