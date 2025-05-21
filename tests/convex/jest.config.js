/**
 * Jest configuration for Convex function tests
 * 
 * This file is only used for running tests and is not loaded
 * by the Convex runtime.
 */

module.exports = {
  preset: 'ts-jest',
  testEnvironment: 'node',
  testMatch: ['**/__tests__/**/*.ts?(x)', '**/?(*.)+(spec|test).ts?(x)'],
  moduleFileExtensions: ['ts', 'tsx', 'js', 'jsx', 'json', 'node'],
  transform: {
    '^.+\\.tsx?$': ['ts-jest', {
      tsconfig: '../../convex/tsconfig.json',
    }],
  },
  setupFilesAfterEnv: ['./setup.ts'],
  // Use this to create proper module name mappings
  moduleNameMapper: {
    '^../(.*)$': '<rootDir>/../../convex/$1'
  },
  // Exclude these files/patterns from being loaded during tests
  modulePathIgnorePatterns: [
    '_generated',
    'node_modules',
    'convex\.json'
  ],
  // Root directory for tests
  rootDir: __dirname,
};