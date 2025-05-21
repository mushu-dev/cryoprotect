import React from 'react';
import Link from 'next/link';

/**
 * MoleculeCard component to display a single molecule in a card format
 * @param {Object} props - Component props
 * @param {Object} props.molecule - Molecule data
 */
export default function MoleculeCard({ molecule }) {
  if (!molecule) {
    return null;
  }

  return (
    <div className="molecule-card" tabIndex="0">
      <h3>{molecule.name || 'Unnamed Molecule'}</h3>
      
      {molecule.pubchem_cid && (
        <div className="molecule-id">PubChem CID: {molecule.pubchem_cid}</div>
      )}
      
      {molecule.formula && (
        <div className="molecule-formula">{molecule.formula}</div>
      )}
      
      <div className="molecule-properties">
        {molecule.molecular_weight && (
          <div className="molecule-property">
            <span className="molecule-property-label">Molecular Weight:</span>
            <span>{molecule.molecular_weight}</span>
          </div>
        )}
        
        {molecule.is_cryoprotectant && (
          <div className="molecule-property">
            <span className="molecule-property-label">Cryoprotectant:</span>
            <span>{molecule.is_cryoprotectant ? 'Yes' : 'No'}</span>
          </div>
        )}
      </div>
      
      <div className="view-details">
        <Link href={`/molecules/${molecule.id}`} className="view-details-link">
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