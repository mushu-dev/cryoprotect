module.exports = {
  // Use jsdom environment for React component tests
  testEnvironment: 'jsdom',
  
  // Extend timeout for async tests
  testTimeout: 20000,
  
  // Test all unit test files
  testMatch: [
    '**/__tests__/**/*.js',
    '**/*.test.js'
  ],
  
  // Ignore test files in node_modules
  testPathIgnorePatterns: [
    '/node_modules/',
    '/.next/',
    '/dist/'
  ],
  
  // Transform files with babel
  transform: {
    '^.+\\.(js|jsx|ts|tsx)$': ['babel-jest', { presets: ['next/babel'] }]
  },
  
  // Handle module aliases - sync with tsconfig paths
  moduleNameMapper: {
    '^@/components/(.*)$': '<rootDir>/src/components/$1',
    '^@/lib/(.*)$': '<rootDir>/src/lib/$1',
    '^@/features/(.*)$': '<rootDir>/src/features/$1',
    '^@/styles/(.*)$': '<rootDir>/src/styles/$1',
    '^@/hooks/(.*)$': '<rootDir>/src/hooks/$1',
    '^@/contexts/(.*)$': '<rootDir>/src/contexts/$1',
    '\\.(css|less|sass|scss)$': '<rootDir>/tests/mocks/styleMock.js',
    '\\.(gif|ttf|eot|svg|png|jpg|jpeg)$': '<rootDir>/tests/mocks/fileMock.js'
  },
  
  // Set up test environment
  setupFilesAfterEnv: [
    '<rootDir>/jest.setup.js'
  ],
  
  // Verbose output for debugging
  verbose: true,
  
  // Enable coverage reporting
  collectCoverage: true,
  collectCoverageFrom: [
    'src/**/*.{js,jsx,ts,tsx}',
    '!src/**/*.d.ts',
    '!src/**/_*.{js,jsx,ts,tsx}',
    '!**/*.stories.{js,jsx,ts,tsx}',
    '!**/node_modules/**'
  ],
  
  // React 18 settings
  testEnvironmentOptions: {
    customExportConditions: [''],
  }
};