/**
 * Unified data fetching hook for CryoProtect.
 * 
 * This hook provides a consistent interface for fetching data from either
 * Convex (when enabled) or the API backend (when Convex is disabled).
 */
import { useQuery } from "convex/react";
import { useState, useEffect } from "react";
import { isEnabled } from "../convex/client";

/**
 * Hook to fetch data from either Convex or API backend.
 * 
 * @param {string} collection - The collection/table name to fetch from
 * @param {object} options - Options for the query
 * @param {string} options.id - Optional ID to fetch a specific item
 * @param {object} options.filters - Optional filters for the query
 * @param {number} options.limit - Optional limit for the query
 * @param {function} options.convexQuery - Function to generate Convex query
 * @returns {object} Data, loading state, and error
 */
export function useData(collection, options = {}) {
  const { id, filters, limit, convexQuery } = options;
  
  // State for API data
  const [apiData, setApiData] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  
  // Try to use Convex if enabled and a convexQuery is provided
  const convexData = isEnabled() && convexQuery 
    ? useQuery(convexQuery, { id, filters, limit }) 
    : undefined;
  
  // Fall back to API if Convex is not enabled or not available
  useEffect(() => {
    if (!isEnabled()) {
      setIsLoading(true);
      
      // Create API endpoint based on collection
      const endpoint = id 
        ? `/api/${collection}/${id}` 
        : `/api/${collection}`;
      
      // Add query parameters if needed
      const queryParams = [];
      if (limit) queryParams.push(`limit=${limit}`);
      if (filters) {
        Object.entries(filters).forEach(([key, value]) => {
          queryParams.push(`${key}=${encodeURIComponent(value)}`);
        });
      }
      
      const url = queryParams.length > 0 
        ? `${endpoint}?${queryParams.join('&')}` 
        : endpoint;
      
      // Fetch from API
      fetch(url)
        .then(res => {
          if (!res.ok) {
            throw new Error(`API returned ${res.status}`);
          }
          return res.json();
        })
        .then(data => {
          setApiData(data);
          setIsLoading(false);
        })
        .catch(err => {
          console.error(`Error fetching from ${url}:`, err);
          setError(err);
          setIsLoading(false);
        });
    }
  }, [collection, id, JSON.stringify(filters), limit]);
  
  return {
    data: isEnabled() ? convexData : apiData,
    isLoading: isEnabled() 
      ? convexData === undefined && !error 
      : isLoading,
    error: isEnabled() ? undefined : error
  };
}

/**
 * Hook to create data in either Convex or API backend.
 * 
 * @param {string} collection - The collection/table name to create in
 * @param {function} options.convexMutation - Function to generate Convex mutation
 * @returns {object} Create function, loading state, and error
 */
export function useCreateData(collection, options = {}) {
  const { convexMutation } = options;
  
  // State for API operations
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  
  // Function to create data
  const createData = async (data) => {
    if (isEnabled() && convexMutation) {
      // Use Convex mutation
      try {
        return await convexMutation(data);
      } catch (err) {
        console.error("Convex mutation error:", err);
        throw err;
      }
    } else {
      // Use API endpoint
      setIsLoading(true);
      setError(null);
      
      try {
        const response = await fetch(`/api/${collection}`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify(data)
        });
        
        if (!response.ok) {
          throw new Error(`API returned ${response.status}`);
        }
        
        const result = await response.json();
        setIsLoading(false);
        return result;
      } catch (err) {
        console.error(`Error creating data in ${collection}:`, err);
        setError(err);
        setIsLoading(false);
        throw err;
      }
    }
  };
  
  return {
    createData,
    isLoading,
    error
  };
}