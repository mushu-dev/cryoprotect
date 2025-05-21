import React from 'react';

/**
 * A component to display molecule properties with visualizations
 * @param {Object} props - Component props
 * @param {Array} props.properties - Array of property objects with name, value, and units
 */
export default function PropertyDisplay({ properties }) {
  if (!properties || properties.length === 0) {
    return (
      <div className="property-display-empty">
        <p>No properties available for this molecule.</p>
      </div>
    );
  }

  // Find min and max values for numeric properties to calculate property bars
  const numericProperties = properties.filter(prop => !isNaN(parseFloat(prop.value)));
  const getPropertyBarWidth = (property) => {
    if (isNaN(parseFloat(property.value))) return 0;
    
    // For demonstration, we'll use a relative scale within the available properties
    const value = parseFloat(property.value);
    const absValue = Math.abs(value);
    
    // Find max absolute value among all properties
    const maxAbsValue = Math.max(...numericProperties.map(p => Math.abs(parseFloat(p.value))));
    
    // Calculate width as percentage (max 95%)
    return (absValue / maxAbsValue) * 95;
  };

  return (
    <div className="property-display">
      <div className="property-list">
        {properties.map((property, index) => {
          const barWidth = getPropertyBarWidth(property);
          const value = parseFloat(property.value);
          const isNegative = value < 0;
          
          return (
            <div className="property-item" key={index}>
              <div className="property-header">
                <div className="property-name">{property.name || 'Property'}</div>
                <div className="property-value">
                  {property.value}
                  {property.units && <span className="property-units">{property.units}</span>}
                </div>
              </div>
              
              {!isNaN(value) && (
                <div className="property-bar-container">
                  <div 
                    className={`property-bar ${isNegative ? 'negative' : 'positive'}`}
                    style={{ width: `${barWidth}%` }}
                  />
                </div>
              )}
            </div>
          );
        })}
      </div>
      
      <style jsx>{`
        .property-display {
          width: 100%;
        }
        
        .property-display-empty {
          padding: 20px;
          text-align: center;
          color: #718096;
        }
        
        .property-list {
          display: flex;
          flex-direction: column;
          gap: 15px;
        }
        
        .property-item {
          width: 100%;
        }
        
        .property-header {
          display: flex;
          justify-content: space-between;
          margin-bottom: 5px;
        }
        
        .property-name {
          font-weight: 500;
          color: #4a5568;
        }
        
        .property-value {
          font-weight: 500;
          color: #2d3748;
        }
        
        .property-units {
          margin-left: 3px;
          color: #718096;
          font-size: 0.9em;
        }
        
        .property-bar-container {
          width: 100%;
          height: 6px;
          background-color: #e2e8f0;
          border-radius: 3px;
          overflow: hidden;
        }
        
        .property-bar {
          height: 100%;
          border-radius: 3px;
          transition: width 0.5s ease;
        }
        
        .property-bar.positive {
          background-color: #4299e1;
        }
        
        .property-bar.negative {
          background-color: #f56565;
        }
      `}</style>
    </div>
  );
}