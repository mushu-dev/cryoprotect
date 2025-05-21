'use server';

import { MOCK_PROTOCOLS } from '../data/mock-protocols';

/**
 * Get a single protocol by ID
 */
export async function getProtocolById(id: string) {
  // In a real application, this would make an API request
  // For now, we'll use the mock data
  const protocol = MOCK_PROTOCOLS.find(p => p.id === id);
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 300));
  
  return protocol;
}

/**
 * Get all protocols
 */
export async function getAllProtocols() {
  // In a real application, this would make an API request
  // For now, we'll use the mock data
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 500));
  
  return MOCK_PROTOCOLS;
}

/**
 * Create a new protocol
 */
export async function createProtocol(data: any) {
  // In a real application, this would make an API request
  console.log('Creating protocol with data:', data);
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 1000));
  
  // Return a mock new protocol
  return {
    id: String(Math.floor(Math.random() * 10000)),
    ...data,
    version: '1.0.0',
    created_at: new Date().toISOString(),
    updated_at: new Date().toISOString()
  };
}

/**
 * Get protocol versions
 */
export async function getProtocolVersions(protocolId: string) {
  // In a real application, this would make an API request
  // For now, we'll return mock data
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 500));
  
  // Generate mock versions
  return [
    {
      id: protocolId,
      version: '1.2.0',
      created_at: '2023-10-15',
      created_by: 'Dr. John Doe',
      description: 'Current version'
    },
    {
      id: `${protocolId}-v1`, 
      version: '1.1.0',
      created_at: '2023-08-15',
      created_by: 'Dr. John Doe',
      description: 'Added improved cooling rates'
    },
    {
      id: `${protocolId}-v0`,
      version: '1.0.0',
      created_at: '2023-07-01',
      created_by: 'Dr. John Doe',
      description: 'Initial version'
    }
  ];
}