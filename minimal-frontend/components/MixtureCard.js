import React from 'react';
import Link from 'next/link';

/**
 * MixtureCard component to display a single mixture in a card format
 * @param {Object} props - Component props
 * @param {Object} props.mixture - Mixture data
 */
export default function MixtureCard({ mixture }) {
  if (!mixture) {
    return null;
  }

  return (
    <div className="mixture-card" tabIndex="0">
      <h3>{mixture.name || 'Unnamed Mixture'}</h3>
      
      {mixture.id && (
        <div className="mixture-id">Mixture ID: {mixture.id}</div>
      )}
      
      {mixture.description && (
        <div className="mixture-description">{mixture.description}</div>
      )}
      
      <div className="mixture-properties">
        {mixture.components && mixture.components.length > 0 && (
          <div className="mixture-property">
            <span className="mixture-property-label">Components:</span>
            <span>{mixture.components.length}</span>
          </div>
        )}
        
        {mixture.freezing_point && (
          <div className="mixture-property">
            <span className="mixture-property-label">Freezing Point:</span>
            <span>{mixture.freezing_point} °C</span>
          </div>
        )}
      </div>
      
      <div className="view-details">
        <Link href={`/mixtures/${mixture.id}`} className="view-details-link">
          View Details
          <svg className="arrow-icon" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
            <path d="M5 12h14"></path>
            <path d="M12 5l7 7-7 7"></path>
          </svg>
        </Link>
      </div>
      
      <style jsx>{`
        .view-details-link {
          display: inline-flex;
          align-items: center;
          gap: 4px;
          padding: 6px 8px;
          transition: all 0.2s ease-in-out;
          border-radius: 4px;
        }
        
        .view-details-link:hover {
          background-color: #f0f7ff;
          text-decoration: none;
          transform: translateX(3px);
        }
        
        .arrow-icon {
          transition: transform 0.2s ease-in-out;
        }
        
        .view-details-link:hover .arrow-icon {
          transform: translateX(2px);
        }
        
        @media (max-width: 480px) {
          .view-details-link {
            width: 100%;
            justify-content: center;
            padding: 8px;
            background-color: #f0f7ff;
            margin-top: 10px;
          }
        }
      `}</style>
    </div>
  );
}