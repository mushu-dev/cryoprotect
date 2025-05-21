'use server';

import { MOCK_EXPERIMENTS } from '../data/mock-experiments';

/**
 * Get a single experiment by ID
 */
export async function getExperimentById(id: string) {
  // In a real application, this would make an API request
  // For now, we'll use the mock data
  const experiment = MOCK_EXPERIMENTS.find(exp => exp.id === id);
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 300));
  
  return experiment;
}

/**
 * Get all experiments
 */
export async function getAllExperiments() {
  // In a real application, this would make an API request
  // For now, we'll use the mock data
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 500));
  
  return MOCK_EXPERIMENTS;
}

/**
 * Create a new experiment
 */
export async function createExperiment(data: any) {
  // In a real application, this would make an API request
  console.log('Creating experiment with data:', data);
  
  // Add a slight delay to simulate network latency
  await new Promise(resolve => setTimeout(resolve, 1000));
  
  // Return a mock new experiment
  return {
    id: String(Math.floor(Math.random() * 10000)),
    ...data,
    date: new Date().toLocaleDateString('en-US', {
      year: 'numeric',
      month: 'long',
      day: 'numeric'
    })
  };
}