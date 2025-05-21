/**
 * Tests for authentication functions
 *
 * Note: These tests can be run without dependency on the client-side 
 * React components as they focus on the server-side logic.
 */

import { createMockContext } from './setup';
import { getCurrentUser, upsertUser, updateUserRole } from '../../convex/auth/users';
import { verifyClerkJWT } from '../../convex/auth/clerk';

// Mark the tests as to be skipped if needed
// Remove the '.skip' when running the tests with all dependencies installed
describe('Authentication Functions', () => {
  // Setup mock context and data
  const ctx = createMockContext();
  const mockUserId = '01234567890123456789012345678901';
  const mockUser = {
    _id: { id: mockUserId, tableName: 'users' },
    _creationTime: Date.now(),
    email: 'test@example.com',
    name: 'Test User',
    role: 'viewer',
    lastLogin: Date.now(),
    createdAt: Date.now(),
    updatedAt: Date.now(),
  };

  beforeEach(() => {
    // Reset mock functions
    jest.resetAllMocks();
    
    // Setup default mocks
    ctx.db.get.mockResolvedValue(mockUser);
    ctx.db.insert.mockResolvedValue({ id: mockUserId, tableName: 'users' });
    ctx.db.patch.mockResolvedValue({ ...mockUser, role: 'scientist' });
    
    // Mock query results
    ctx.db.query().filter().first.mockResolvedValue(mockUser);
    ctx.db.query().filter().collect.mockResolvedValue([mockUser]);
    ctx.db.query().collect.mockResolvedValue([mockUser]);
  });

  describe('getCurrentUser', () => {
    test('returns authenticated user when logged in', async () => {
      const result = await getCurrentUser.handler(ctx as any, {});
      
      expect(result.isAuthenticated).toBe(true);
      expect(result.user).toEqual(mockUser);
      expect(result.roles).toEqual(['viewer']);
    });

    test('returns null user when not authenticated', async () => {
      // Mock unauthenticated context
      const unauthCtx = {
        ...ctx,
        auth: {
          userId: null,
          subject: null,
          isAuthenticated: false,
        },
      };
      
      const result = await getCurrentUser.handler(unauthCtx as any, {});
      
      expect(result.isAuthenticated).toBe(false);
      expect(result.user).toBeNull();
      expect(result.roles).toEqual([]);
    });
  });

  describe('upsertUser', () => {
    test('creates a new user if not found', async () => {
      // Mock that user doesn't exist
      ctx.db.query().filter().first.mockResolvedValue(null);
      
      const args = {
        email: 'test@example.com',
        name: 'New User',
      };
      
      const result = await upsertUser.handler(ctx as any, args);
      
      expect(ctx.db.insert).toHaveBeenCalledWith('users', expect.objectContaining({
        email: 'test@example.com',
        name: 'New User',
        role: 'viewer',
      }));
      expect(result).toEqual({ id: mockUserId, tableName: 'users' });
    });

    test('updates existing user if found', async () => {
      // Mock that user exists
      ctx.db.query().filter().first.mockResolvedValue(mockUser);
      
      const args = {
        email: 'test@example.com',
        name: 'Updated User',
      };
      
      const result = await upsertUser.handler(ctx as any, args);
      
      expect(ctx.db.patch).toHaveBeenCalledWith(
        mockUser._id,
        expect.objectContaining({
          name: 'Updated User',
        })
      );
      expect(result).toEqual(mockUser._id);
    });
  });

  describe('updateUserRole', () => {
    test('updates user role as admin', async () => {
      // Mock that current user is admin
      ctx.db.query().filter().first.mockResolvedValue({
        ...mockUser,
        role: 'admin',
      });
      
      const args = {
        userId: { id: mockUserId, tableName: 'users' },
        role: 'scientist',
      };
      
      const result = await updateUserRole.handler(ctx as any, args);
      
      expect(ctx.db.patch).toHaveBeenCalledWith(
        { id: mockUserId, tableName: 'users' },
        expect.objectContaining({
          role: 'scientist',
        })
      );
      expect(result.role).toBe('scientist');
    });
  });

  describe('verifyClerkJWT', () => {
    test('returns valid token information', async () => {
      const args = {
        clerkJWT: 'valid_jwt_token',
      };
      
      const result = await verifyClerkJWT.handler(ctx as any, args);
      
      expect(result.isValid).toBe(true);
      expect(result.subject).toBe('user@example.com');
    });
  });
});