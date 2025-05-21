import { httpRouter } from "convex/server";
import { httpAction } from "../_generated/server";
import { v } from "convex/values";

/**
 * HTTP router for RESTful API endpoints that match the Supabase-like
 * interface expected by the ConvexAdapter in Python.
 */
const http = httpRouter();

/**
 * Query endpoint - handles database queries
 * Expected format: { "table": string, "select": string[] | "*", "where": {}, "limit": number }
 */
http.route({
  path: "/api/query",
  method: "POST",
  handler: httpAction({
    handler: async (ctx, request) => {
      const body = await request.json();
      const { table, select, where, limit } = body;
      
      let query = ctx.db.query(table);
      
      // Apply where conditions if provided
      if (where) {
        Object.entries(where).forEach(([key, value]) => {
          if (value !== null && value !== undefined) {
            query = query.filter(q => q.eq(key, value));
          }
        });
      }
      
      // Apply limit if provided
      if (limit) {
        query = query.take(limit);
      } else {
        query = query.take(1000); // Default limit
      }
      
      const results = await query.collect();
      
      return new Response(JSON.stringify({ data: results }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    },
  }),
});

/**
 * Insert endpoint - handles data insertion
 * Expected format: { "table": string, "data": object | object[] }
 */
http.route({
  path: "/api/insert",
  method: "POST",
  handler: httpAction({
    handler: async (ctx, request) => {
      const body = await request.json();
      const { table, data } = body;
      
      let insertedIds = [];
      
      if (Array.isArray(data)) {
        // Batch insert
        for (const item of data) {
          const id = await ctx.db.insert(table, item);
          insertedIds.push(id);
        }
      } else {
        // Single insert
        const id = await ctx.db.insert(table, data);
        insertedIds.push(id);
      }
      
      return new Response(JSON.stringify({ 
        data: insertedIds,
        success: true,
        status: 201
      }), {
        status: 201,
        headers: { "Content-Type": "application/json" }
      });
    },
  }),
});

/**
 * Update endpoint - handles data updates
 * Expected format: { "table": string, "id": string, "data": object }
 */
http.route({
  path: "/api/update",
  method: "POST",
  handler: httpAction({
    handler: async (ctx, request) => {
      const body = await request.json();
      const { table, id, data } = body;
      
      await ctx.db.patch(id, data);
      
      return new Response(JSON.stringify({ 
        success: true,
        status: 200
      }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    },
  }),
});

/**
 * Delete endpoint - handles data deletion
 * Expected format: { "table": string, "id": string }
 */
http.route({
  path: "/api/delete",
  method: "POST",
  handler: httpAction({
    handler: async (ctx, request) => {
      const body = await request.json();
      const { table, id } = body;
      
      await ctx.db.delete(id);
      
      return new Response(JSON.stringify({ 
        success: true,
        status: 200
      }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    },
  }),
});

/**
 * Mock auth endpoint - provides authentication compatibility
 */
http.route({
  path: "/api/auth",
  method: "POST",
  handler: httpAction({
    handler: async (ctx, request) => {
      // This is a mock endpoint that returns success for compatibility
      return new Response(JSON.stringify({ 
        user: { id: "convex-user-1" },
        session: { access_token: "mock-token" },
        success: true,
        status: 200
      }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    },
  }),
});

/**
 * Health check endpoint
 */
http.route({
  path: "/api/health",
  method: "GET",
  handler: httpAction({
    handler: async (ctx, request) => {
      return new Response(JSON.stringify({ 
        status: "ok",
        message: "Convex HTTP API is running"
      }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    },
  }),
});

export default http;