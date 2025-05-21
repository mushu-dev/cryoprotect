/**
 * Test setup for Convex function tests
 * 
 * NOTE: This file is only used for jest testing and is not meant
 * to be loaded by the Convex runtime.
 */

// This conditional check prevents this file from being executed
// in a Convex environment when pushing deployment config
if (typeof jest !== 'undefined') {
  // Mock environment variables
  process.env.CLERK_PUBLISHABLE_KEY = 'test_clerk_key';
  process.env.CONVEX_URL = 'https://test.convex.cloud';

  // Set timeout for tests
  jest.setTimeout(10000);
}

// Safely define ConvexError if it doesn't exist
if (typeof global.ConvexError === 'undefined') {
  global.ConvexError = class ConvexError extends Error {
    constructor(public readonly code: string, message: string) {
      super(message);
      this.name = 'ConvexError';
    }
  };
}

/**
 * Helper functions for tests
 * This is only used in test environments and creates mock objects
 * that simulate the Convex context.
 */
export const createMockContext = () => {
  // Only create mock functions if jest is defined
  const createMockFn = typeof jest !== 'undefined' 
    ? () => jest.fn()
    : () => (...args: any[]) => undefined;
    
  return {
    auth: {
      userId: 'test_user_id',
      subject: 'test@example.com',
      isAuthenticated: true,
    },
    db: {
      get: createMockFn(),
      query: createMockFn().mockReturnValue({
        filter: createMockFn().mockReturnValue({
          collect: createMockFn(),
          first: createMockFn(),
        }),
        withIndex: createMockFn().mockReturnValue({
          filter: createMockFn().mockReturnValue({
            collect: createMockFn(),
            first: createMockFn(),
          }),
          eq: createMockFn().mockReturnValue({
            collect: createMockFn(),
            first: createMockFn(),
          }),
        }),
        collect: createMockFn(),
        first: createMockFn(),
      }),
      insert: createMockFn(),
      patch: createMockFn(),
      delete: createMockFn(),
    },
  };
};