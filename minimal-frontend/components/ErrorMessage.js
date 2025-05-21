import React from 'react';

/**
 * Error message component to display when API requests fail
 * @param {Object} props - Component props
 * @param {string} props.message - Error message to display
 * @param {Function} props.onRetry - Function to call when retry button is clicked
 */
export default function ErrorMessage({ message, onRetry }) {
  return (
    <div className="error-message">
      <h3>Error</h3>
      <p>{message || 'An error occurred while fetching data.'}</p>
      {onRetry && (
        <button className="button" onClick={onRetry}>
          Retry
        </button>
      )}
    </div>
  );
}