import { useState, useEffect } from 'react';
import { getProtocolService } from '../../../features/common/api/service-factory';

// Basic protocol type
const mockProtocols = [
  {
    id: '1',
    name: 'Standard Cryopreservation Protocol',
    description: 'A standard protocol for cryopreservation of biological samples.',
    category: 'Cryopreservation',
    steps: [
      { title: 'Sample Preparation', description: 'Prepare the sample for cryopreservation' },
      { title: 'Add Cryoprotectant', description: 'Add appropriate cryoprotectant solution' },
      { title: 'Controlled Cooling', description: 'Cool the sample at controlled rate' },
      { title: 'Storage', description: 'Store in appropriate conditions' }
    ],
    createdAt: new Date(Date.now() - 60 * 24 * 60 * 60 * 1000),
    author: 'Dr. Smith'
  },
  {
    id: '2',
    name: 'Vitrification Protocol',
    description: 'Rapid cooling protocol to achieve vitrification of samples.',
    category: 'Vitrification',
    steps: [
      { title: 'Sample Preparation', description: 'Prepare the sample for vitrification' },
      { title: 'Add Cryoprotectant Mixture', description: 'Add vitrification solution' },
      { title: 'Rapid Cooling', description: 'Plunge into liquid nitrogen' },
      { title: 'Storage', description: 'Store in liquid nitrogen' }
    ],
    createdAt: new Date(Date.now() - 45 * 24 * 60 * 60 * 1000),
    author: 'Dr. Johnson'
  }
];

// Now using service factory instead of direct instantiation
export function useProtocols(initialParams = {}) {
  const [protocols, setProtocols] = useState([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [params, setParams] = useState(initialParams)
  
  useEffect(() => {
    const fetchProtocols = async () => {
      try {
        setLoading(true)
        
        // In a real app, this would call an API
        // const service = serviceFactory.getProtocolService()
        // const data = await service.getProtocols(params)
        
        // For now, return mock data
        setTimeout(() => {
          setProtocols(mockProtocols)
          setLoading(false)
        }, 800)
        
      } catch (err) {
        console.error('Error fetching protocols:', err)
        setError(err instanceof Error ? err : new Error('Failed to fetch protocols'))
        setLoading(false)
      }
    }
    
    fetchProtocols()
  }, [params])
  
  // Function to update search parameters
  const updateParams = (newParams) => {
    setParams(prev => ({ ...prev, ...newParams }))
  }
  
  // Function to get a specific protocol by ID
  const getProtocolById = async (id) => {
    try {
      // In a real app, this would call an API
      // const service = serviceFactory.getProtocolService()
      // return await service.getProtocolById(id)
      
      // For now, return mock data
      return mockProtocols.find(p => p.id === id) || null
    } catch (err) {
      console.error(`Error fetching protocol with ID ${id}:`, err)
      throw err instanceof Error ? err : new Error(`Failed to fetch protocol with ID ${id}`)
    }
  }
  
  // Function to create a new protocol
  const createProtocol = async (protocolData) => {
    try {
      // In a real app, this would call an API
      // const service = serviceFactory.getProtocolService()
      // return await service.createProtocol(protocolData)
      
      // For now, just log and return mock data
      console.log('Creating protocol with data:', protocolData)
      const newProtocol = {
        id: String(Date.now()),
        ...protocolData,
        createdAt: new Date(),
        updatedAt: new Date()
      }
      
      // Add to local state
      setProtocols(prev => [newProtocol, ...prev])
      
      return newProtocol
    } catch (err) {
      console.error('Error creating protocol:', err)
      throw err instanceof Error ? err : new Error('Failed to create protocol')
    }
  }
  
  // Function to update an existing protocol
  const updateProtocol = async (id, updates) => {
    try {
      // In a real app, this would call an API
      // const service = serviceFactory.getProtocolService()
      // return await service.updateProtocol(id, updates)
      
      // For now, just update local state
      console.log(`Updating protocol ${id} with:`, updates)
      const updatedProtocols = protocols.map(p => 
        p.id === id ? { ...p, ...updates, updatedAt: new Date() } : p
      )
      
      setProtocols(updatedProtocols)
      
      return updatedProtocols.find(p => p.id === id)
    } catch (err) {
      console.error(`Error updating protocol with ID ${id}:`, err)
      throw err instanceof Error ? err : new Error(`Failed to update protocol with ID ${id}`)
    }
  }
  
  // Function to delete a protocol
  const deleteProtocol = async (id) => {
    try {
      // In a real app, this would call an API
      // const service = serviceFactory.getProtocolService()
      // await service.deleteProtocol(id)
      
      // For now, just update local state
      console.log(`Deleting protocol ${id}`)
      setProtocols(prev => prev.filter(p => p.id !== id))
      
      return true
    } catch (err) {
      console.error(`Error deleting protocol with ID ${id}:`, err)
      throw err instanceof Error ? err : new Error(`Failed to delete protocol with ID ${id}`)
    }
  }
  
  return {
    protocols,
    loading,
    error,
    updateParams,
    getProtocolById,
    createProtocol,
    updateProtocol,
    deleteProtocol
  }
}