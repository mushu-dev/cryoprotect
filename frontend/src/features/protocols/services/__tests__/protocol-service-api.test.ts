/**
 * Tests for the Protocol API Service
 */

import { ApiProtocolService } from '../protocol-service-api';
import { Protocol } from '../protocol-service';

// Mock fetch globally
global.fetch = jest.fn();

describe('ApiProtocolService', () => {
  let service: ApiProtocolService;
  
  beforeEach(() => {
    jest.clearAllMocks();
    service = new ApiProtocolService('/api/v1');
  });
  
  it('should initialize with the correct API URL', () => {
    expect((service as any).apiUrl).toBe('/api/v1/protocols');
  });
  
  it('should fetch protocols with appropriate parameters', async () => {
    // Mock successful response
    const mockResponse = {
      data: [{ id: '123', name: 'Test Protocol' }],
      total: 1,
      page: 1,
      per_page: 10,
      total_pages: 1
    };
    
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: true,
      json: jest.fn().mockResolvedValueOnce(mockResponse)
    });
    
    const result = await service.getProtocols({ page: 1, per_page: 10 });
    
    expect(global.fetch).toHaveBeenCalledWith(
      '/api/v1/protocols?page=1&per_page=10',
      expect.objectContaining({ 
        headers: expect.objectContaining({
          'Content-Type': 'application/json'
        })
      })
    );
    
    expect(result).toEqual(mockResponse);
  });
  
  it('should get a single protocol by ID', async () => {
    const mockProtocol = {
      id: '123',
      name: 'Test Protocol',
      steps: []
    };
    
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: true,
      json: jest.fn().mockResolvedValueOnce(mockProtocol)
    });
    
    const result = await service.getProtocol('123' as any);
    
    expect(global.fetch).toHaveBeenCalledWith(
      '/api/v1/protocols/123',
      expect.objectContaining({ 
        headers: expect.objectContaining({
          'Content-Type': 'application/json'
        })
      })
    );
    
    expect(result).toEqual(mockProtocol);
  });
  
  it('should create a protocol', async () => {
    const newProtocol = {
      name: 'New Protocol',
      steps: [],
      target_concentration: 10,
      sample_type: 'cell_line',
      starting_temperature: 4,
      created_by: 'user123'
    };
    
    const mockResponse = {
      ...newProtocol,
      id: '123',
      created_at: '2023-01-01T00:00:00Z',
      updated_at: '2023-01-01T00:00:00Z'
    };
    
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: true,
      json: jest.fn().mockResolvedValueOnce(mockResponse)
    });
    
    const result = await service.createProtocol(newProtocol as any);
    
    expect(global.fetch).toHaveBeenCalledWith(
      '/api/v1/protocols',
      expect.objectContaining({ 
        method: 'POST',
        headers: expect.objectContaining({
          'Content-Type': 'application/json'
        }),
        body: JSON.stringify(newProtocol)
      })
    );
    
    expect(result).toEqual(mockResponse);
  });
  
  it('should handle API errors gracefully', async () => {
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: false,
      status: 404,
      statusText: 'Not Found',
      text: jest.fn().mockResolvedValueOnce('{"error": "Protocol not found"}')
    });
    
    await expect(service.getProtocol('non-existent' as any))
      .rejects
      .toThrow('Protocol not found');
  });
  
  it('should use authentication when provided', async () => {
    const authenticatedService = new ApiProtocolService('/api/v1', 'test-token');
    
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: true,
      json: jest.fn().mockResolvedValueOnce({})
    });
    
    await authenticatedService.getProtocols();
    
    expect(global.fetch).toHaveBeenCalledWith(
      expect.any(String),
      expect.objectContaining({ 
        headers: expect.objectContaining({
          'Authorization': 'Bearer test-token'
        })
      })
    );
  });
  
  it('should design a protocol for a mixture', async () => {
    const mockParams = {
      target_concentration: 10,
      sample_type: 'cell_line',
      starting_temperature: 4
    };
    
    const mockProtocol = {
      id: '123',
      name: 'Generated Protocol',
      steps: []
    };
    
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: true,
      json: jest.fn().mockResolvedValueOnce(mockProtocol)
    });
    
    const result = await service.designProtocol('mixture-123' as any, mockParams);
    
    expect(global.fetch).toHaveBeenCalledWith(
      '/api/v1/protocols/design/mixture-123',
      expect.objectContaining({ 
        method: 'POST',
        headers: expect.objectContaining({
          'Content-Type': 'application/json'
        }),
        body: JSON.stringify(mockParams)
      })
    );
    
    expect(result).toEqual(mockProtocol);
  });
  
  it('should save a protocol for a mixture', async () => {
    const mockProtocolData = {
      name: 'Saved Protocol',
      target_concentration: 10,
      sample_type: 'cell_line',
      starting_temperature: 4,
      steps: []
    };
    
    const mockSavedProtocol = {
      ...mockProtocolData,
      id: '123',
      created_at: '2023-01-01T00:00:00Z'
    };
    
    (global.fetch as jest.Mock).mockResolvedValueOnce({
      ok: true,
      json: jest.fn().mockResolvedValueOnce(mockSavedProtocol)
    });
    
    const result = await service.saveProtocol('mixture-123' as any, mockProtocolData);
    
    expect(global.fetch).toHaveBeenCalledWith(
      '/api/v1/protocols/save/mixture-123',
      expect.objectContaining({ 
        method: 'POST',
        headers: expect.objectContaining({
          'Content-Type': 'application/json'
        }),
        body: JSON.stringify(mockProtocolData)
      })
    );
    
    expect(result).toEqual(mockSavedProtocol);
  });
});