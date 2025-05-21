import { useState, useEffect } from 'react';
import { ConvexProvider, ConvexReactClient } from 'convex/react';
import { convex } from './client';

/**
 * ConvexClientProvider component
 * 
 * This component wraps the application with the Convex client provider,
 * handling authentication and client lifecycle.
 * 
 * @param {Object} props - Component props
 * @param {React.ReactNode} props.children - Child components
 */
export function ConvexClientProvider({ children }) {
  const [isAuthenticated, setIsAuthenticated] = useState(false);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);
  
  // Set up anonymous authentication when the component mounts
  useEffect(() => {
    const setupAuth = async () => {
      setIsLoading(true);
      
      try {
        // Set up anonymous authentication
        await convex.setAuth({ type: 'anonymous' });
        setIsAuthenticated(true);
        setError(null);
        console.log('Convex client authenticated (anonymous)');
      } catch (error) {
        console.error('Convex authentication error:', error);
        setIsAuthenticated(false);
        setError(error);
      } finally {
        setIsLoading(false);
      }
    };
    
    setupAuth();
    
    // Clean up authentication when the component unmounts
    return () => {
      try {
        convex.clearAuth();
        console.log('Convex authentication cleared');
      } catch (error) {
        console.error('Error clearing Convex auth:', error);
      }
    };
  }, []);
  
  // Show a loading message if authentication is still in progress
  if (isLoading) {
    return (
      <div className="convex-loading">
        <p>Connecting to Convex...</p>
      </div>
    );
  }
  
  // Show an error message if authentication failed
  if (error) {
    return (
      <div className="convex-error">
        <h3>Error connecting to Convex</h3>
        <p>{error.message || 'Unknown error'}</p>
        <style jsx>{`
          .convex-error {
            padding: 15px;
            margin: 20px 0;
            border: 1px solid #f5c6cb;
            border-radius: 4px;
            color: #721c24;
            background-color: #f8d7da;
          }
        `}</style>
      </div>
    );
  }
  
  // Render the children with the Convex provider
  return (
    <ConvexProvider client={convex}>
      {children}
      <style jsx>{`
        .convex-loading {
          display: flex;
          justify-content: center;
          align-items: center;
          height: 100px;
          font-weight: 500;
          color: #3182ce;
        }
      `}</style>
    </ConvexProvider>
  );
}