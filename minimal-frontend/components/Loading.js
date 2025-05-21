import React from 'react';

/**
 * Loading component to display when fetching data
 */
export default function Loading() {
  return (
    <div className="loading-container">
      <div className="loading-spinner"></div>
      <p>Loading data...</p>
    </div>
  );
}