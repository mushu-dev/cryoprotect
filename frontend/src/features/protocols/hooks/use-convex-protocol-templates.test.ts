/**
 * Tests for Protocol Template Management System hooks
 */

import { renderHook, act } from '@testing-library/react-hooks';
import { useProtocolTemplates, useTemplateVersions, useTemplateComparison } from './use-convex-protocol-templates';

// Mock the Convex hooks
jest.mock('convex/react', () => ({
  useQuery: jest.fn(),
  useMutation: jest.fn(() => jest.fn())
}));

// Mock the error handler
jest.mock('../../errors/use-error-handler', () => ({
  useErrorHandler: jest.fn(() => ({
    handleError: jest.fn()
  }))
}));

// Mock the toast
jest.mock('../../ui/use-toast', () => ({
  useToast: jest.fn(() => ({
    toast: jest.fn()
  }))
}));

describe('Protocol Template Hooks', () => {
  describe('useProtocolTemplates', () => {
    it('should return templates and management functions', () => {
      const mockTemplates = [
        { _id: '1', name: 'Template 1', isTemplate: true },
        { _id: '2', name: 'Template 2', isTemplate: true }
      ];
      
      // Mock the useQuery to return our mock templates
      const useQuery = require('convex/react').useQuery;
      useQuery.mockReturnValue(mockTemplates);
      
      const { result } = renderHook(() => useProtocolTemplates());
      
      expect(result.current.templates).toEqual(mockTemplates);
      expect(typeof result.current.createTemplate).toBe('function');
      expect(typeof result.current.updateTemplate).toBe('function');
      expect(typeof result.current.deleteTemplate).toBe('function');
      expect(typeof result.current.createFromTemplate).toBe('function');
      expect(typeof result.current.getTemplateCategories).toBe('function');
    });
    
    it('should return loading state when templates are undefined', () => {
      // Mock the useQuery to return undefined (loading state)
      const useQuery = require('convex/react').useQuery;
      useQuery.mockReturnValue(undefined);
      
      const { result } = renderHook(() => useProtocolTemplates());
      
      expect(result.current.templates).toEqual([]);
      expect(result.current.loading).toBe(true);
    });
    
    it('should provide correct template categories', () => {
      // Mock the useQuery
      const useQuery = require('convex/react').useQuery;
      useQuery.mockReturnValue([]);
      
      const { result } = renderHook(() => useProtocolTemplates());
      
      const categories = result.current.getTemplateCategories();
      
      expect(categories).toContain('Vitrification');
      expect(categories).toContain('Slow Freezing');
      expect(categories).toContain('Sample Preparation');
    });
  });
  
  describe('useTemplateVersions', () => {
    it('should return template versions when templateName is provided', () => {
      const mockVersions = [
        { _id: '1', name: 'Template 1', version: 1 },
        { _id: '2', name: 'Template 1', version: 2 }
      ];
      
      // Mock the useQuery to return our mock versions
      const useQuery = require('convex/react').useQuery;
      useQuery.mockReturnValue(mockVersions);
      
      const { result } = renderHook(() => useTemplateVersions('Template 1'));
      
      expect(result.current.versions).toEqual(mockVersions);
      expect(result.current.loading).toBe(false);
    });
    
    it('should not query when templateName is not provided', () => {
      // Mock the useQuery
      const useQuery = require('convex/react').useQuery;
      
      renderHook(() => useTemplateVersions());
      
      // Expect the query to have been called with "skip"
      expect(useQuery).toHaveBeenCalledWith(expect.anything(), "skip");
    });
  });
  
  describe('useTemplateComparison', () => {
    it('should provide compareTemplates function', () => {
      global.fetch = jest.fn().mockImplementation(() => 
        Promise.resolve({
          ok: true,
          json: () => Promise.resolve({})
        })
      );
      
      const { result } = renderHook(() => useTemplateComparison());
      
      expect(typeof result.current.compareTemplates).toBe('function');
    });
  });
});