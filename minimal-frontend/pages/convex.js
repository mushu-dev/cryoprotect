import React, { useState } from 'react';
import Layout from '../components/Layout';
import { ConvexMoleculesList } from '../src/components/convex/ConvexMoleculesList';
import { ConvexMixturesList } from '../src/components/convex/ConvexMixturesList';
import { 
  useMolecules, 
  useAddMolecule, 
  useUpdateMolecule, 
  useDeleteMolecule 
} from '../src/convex/hooks';

export default function ConvexDemoPage() {
  const [activeTab, setActiveTab] = useState('molecules');
  const [newMolecule, setNewMolecule] = useState({
    name: '',
    formula: '',
    pubchemCid: '',
    description: ''
  });
  
  // Get molecules
  const molecules = useMolecules();
  
  // Get mutations
  const addMolecule = useAddMolecule();
  const updateMolecule = useUpdateMolecule();
  const deleteMolecule = useDeleteMolecule();
  
  // Handler for adding a new molecule
  const handleAddMolecule = async (e) => {
    e.preventDefault();
    try {
      await addMolecule(newMolecule);
      // Clear form after successful addition
      setNewMolecule({
        name: '',
        formula: '',
        pubchemCid: '',
        description: ''
      });
    } catch (error) {
      console.error('Error adding molecule:', error);
    }
  };
  
  // Handler for form field changes
  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setNewMolecule(prev => ({
      ...prev,
      [name]: value
    }));
  };
  
  return (
    <Layout title="Convex Integration">
      <div className="page-header">
        <h1 className="title">Convex Integration Demo</h1>
        <p className="description">
          This page demonstrates the integration with Convex, our real-time database.
        </p>
      </div>
      
      <div className="card">
        <div className="tabs">
          <button 
            className={`tab ${activeTab === 'molecules' ? 'active' : ''}`}
            onClick={() => setActiveTab('molecules')}
          >
            Molecules
          </button>
          <button 
            className={`tab ${activeTab === 'mixtures' ? 'active' : ''}`}
            onClick={() => setActiveTab('mixtures')}
          >
            Mixtures
          </button>
          <button 
            className={`tab ${activeTab === 'add' ? 'active' : ''}`}
            onClick={() => setActiveTab('add')}
          >
            Add New
          </button>
        </div>
        
        <div className="tab-content">
          {activeTab === 'molecules' && (
            <div>
              <h2>Molecules from Convex</h2>
              <p>These molecules are fetched directly from the Convex real-time database.</p>
              <ConvexMoleculesList />
            </div>
          )}
          
          {activeTab === 'mixtures' && (
            <div>
              <h2>Mixtures from Convex</h2>
              <p>These mixtures are fetched directly from the Convex real-time database.</p>
              <ConvexMixturesList />
            </div>
          )}
          
          {activeTab === 'add' && (
            <div>
              <h2>Add New Molecule</h2>
              <p>Add a new molecule to the Convex database.</p>
              
              <form onSubmit={handleAddMolecule} className="add-form">
                <div className="form-group">
                  <label htmlFor="name">Name:</label>
                  <input
                    type="text"
                    id="name"
                    name="name"
                    value={newMolecule.name}
                    onChange={handleInputChange}
                    required
                  />
                </div>
                
                <div className="form-group">
                  <label htmlFor="formula">Formula:</label>
                  <input
                    type="text"
                    id="formula"
                    name="formula"
                    value={newMolecule.formula}
                    onChange={handleInputChange}
                  />
                </div>
                
                <div className="form-group">
                  <label htmlFor="pubchemCid">PubChem CID:</label>
                  <input
                    type="text"
                    id="pubchemCid"
                    name="pubchemCid"
                    value={newMolecule.pubchemCid}
                    onChange={handleInputChange}
                  />
                </div>
                
                <div className="form-group">
                  <label htmlFor="description">Description:</label>
                  <textarea
                    id="description"
                    name="description"
                    value={newMolecule.description}
                    onChange={handleInputChange}
                    rows="4"
                  />
                </div>
                
                <button type="submit" className="button">Add Molecule</button>
              </form>
            </div>
          )}
        </div>
      </div>
      
      <div className="card" style={{ marginTop: '30px' }}>
        <h2 className="subtitle">Convex Features</h2>
        <div className="features-grid">
          <div className="feature-item">
            <h3>Real-time Updates</h3>
            <p>
              Changes to the data are automatically synchronized across all clients.
              Try opening this page in multiple windows to see changes in real-time.
            </p>
          </div>
          
          <div className="feature-item">
            <h3>Optimistic Updates</h3>
            <p>
              The UI updates immediately when you perform actions, then syncs with the server.
              This provides a responsive user experience even with network latency.
            </p>
          </div>
          
          <div className="feature-item">
            <h3>Automatic Caching</h3>
            <p>
              Convex automatically caches query results and intelligently invalidates
              the cache when data changes, reducing unnecessary network requests.
            </p>
          </div>
          
          <div className="feature-item">
            <h3>Type Safety</h3>
            <p>
              Our integration with Convex includes TypeScript types for all queries
              and mutations, ensuring type safety throughout the application.
            </p>
          </div>
        </div>
      </div>
      
      <style jsx>{`
        .tabs {
          display: flex;
          border-bottom: 1px solid #ddd;
          margin-bottom: 20px;
        }
        
        .tab {
          padding: 10px 20px;
          background: none;
          border: none;
          border-bottom: 2px solid transparent;
          cursor: pointer;
          font-size: 16px;
          font-weight: 500;
        }
        
        .tab.active {
          border-bottom: 2px solid #3182ce;
          color: #3182ce;
        }
        
        .tab-content {
          padding: 10px;
        }
        
        .add-form {
          max-width: 600px;
          margin: 0 auto;
        }
        
        .form-group {
          margin-bottom: 15px;
        }
        
        .form-group label {
          display: block;
          margin-bottom: 5px;
          font-weight: 500;
        }
        
        .form-group input,
        .form-group textarea {
          width: 100%;
          padding: 8px;
          border: 1px solid #ddd;
          border-radius: 4px;
        }
        
        .button {
          background-color: #3182ce;
          color: white;
          border: none;
          padding: 10px 15px;
          border-radius: 4px;
          cursor: pointer;
          font-weight: 500;
        }
        
        .button:hover {
          background-color: #2c5282;
        }
        
        .features-grid {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
          gap: 20px;
          margin-top: 20px;
        }
        
        .feature-item {
          padding: 15px;
          border: 1px solid #eee;
          border-radius: 5px;
          background-color: #f9f9f9;
        }
        
        .feature-item h3 {
          margin-top: 0;
          color: #3182ce;
        }
      `}</style>
    </Layout>
  );
}