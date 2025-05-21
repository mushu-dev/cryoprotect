import { httpRouter } from "convex/server";
import { httpAction } from "../_generated/server";
import { v } from "convex/values";

/**
 * HTTP router for our RESTful API
 */
const http = httpRouter();

/**
 * Query endpoint
 */
http.route({
  path: "/api/query",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      const body = await request.json();
      const { table, filters, columns, order, limit, offset } = body;
      
      // Start building the query
      let query = ctx.db.query(table);
      
      // Apply filters if provided
      if (filters && Object.keys(filters).length > 0) {
        query = query.filter(q => {
          // Start with a condition that will always be true
          let conditions = q.eq(1, 1);
          
          // Add conditions for each filter
          for (const [field, value] of Object.entries(filters)) {
            conditions = q.and(conditions, q.eq(q.field(field), value));
          }
          
          return conditions;
        });
      }
      
      // Execute the query to get the results
      let results = await query.collect();
      
      // Apply pagination
      if (offset && offset > 0) {
        results = results.slice(offset);
      }
      
      if (limit && limit > 0) {
        results = results.slice(0, limit);
      }
      
      // Filter columns if specified
      if (columns && columns !== '*') {
        const columnsToInclude = columns.split(',').map(c => c.trim());
        results = results.map(record => {
          const filteredRecord: Record<string, any> = {};
          columnsToInclude.forEach(col => {
            if (record.hasOwnProperty(col)) {
              filteredRecord[col] = record[col];
            }
          });
          return filteredRecord;
        });
      }
      
      return new Response(JSON.stringify({ data: results, error: null }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error executing query:', error);
      return new Response(JSON.stringify({ 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  }),
});

/**
 * Insert endpoint
 */
http.route({
  path: "/api/insert",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      const body = await request.json();
      const { table, data } = body;
      
      // Add timestamps
      const now = Date.now();
      
      // Handle both single object and array of objects
      let dataToInsert = Array.isArray(data) ? data : [data];
      
      // Add timestamps to each record
      dataToInsert = dataToInsert.map(record => ({
        ...record,
        createdAt: now,
        updatedAt: now
      }));
      
      // Insert the data
      const ids = await Promise.all(
        dataToInsert.map(record => ctx.db.insert(table, record))
      );
      
      // Return the inserted records with their IDs
      const results = dataToInsert.map((record, index) => ({
        ...record,
        id: ids[index]
      }));
      
      return new Response(JSON.stringify({ data: results, error: null }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error inserting data:', error);
      return new Response(JSON.stringify({ 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  })
});

/**
 * Update endpoint
 */
http.route({
  path: "/api/update",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      const body = await request.json();
      const { table, data, filters } = body;
      
      // Add updated timestamp
      const updateData = {
        ...data,
        updatedAt: Date.now()
      };
      
      // Find records matching the filters
      let query = ctx.db.query(table);
      
      if (filters && Object.keys(filters).length > 0) {
        query = query.filter(q => {
          let conditions = q.eq(1, 1);
          
          for (const [field, value] of Object.entries(filters)) {
            conditions = q.and(conditions, q.eq(q.field(field), value));
          }
          
          return conditions;
        });
      }
      
      const records = await query.collect();
      
      // Update each record
      const updatedRecords = await Promise.all(
        records.map(async record => {
          const id = record._id;
          await ctx.db.patch(table, id, updateData);
          return { ...record, ...updateData, id };
        })
      );
      
      return new Response(JSON.stringify({ data: updatedRecords, error: null }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error updating data:', error);
      return new Response(JSON.stringify({ 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  })
});

/**
 * Delete endpoint
 */
http.route({
  path: "/api/delete",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      const body = await request.json();
      const { table, filters } = body;
      
      // Find records matching the filters
      let query = ctx.db.query(table);
      
      if (filters && Object.keys(filters).length > 0) {
        query = query.filter(q => {
          let conditions = q.eq(1, 1);
          
          for (const [field, value] of Object.entries(filters)) {
            conditions = q.and(conditions, q.eq(q.field(field), value));
          }
          
          return conditions;
        });
      }
      
      const records = await query.collect();
      
      // Delete each record
      await Promise.all(
        records.map(record => ctx.db.delete(table, record._id))
      );
      
      return new Response(JSON.stringify({ data: records, error: null }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error deleting data:', error);
      return new Response(JSON.stringify({ 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  })
});

/**
 * Auth signin endpoint
 */
http.route({
  path: "/api/auth/signin",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      const body = await request.json();
      const { email, password } = body;
      
      // Find the user
      const users = await ctx.db
        .query("users")
        .filter(q => q.eq(q.field("email"), email))
        .collect();
      
      if (users.length === 0) {
        return new Response(JSON.stringify({ 
          data: null, 
          error: "Invalid email or password"
        }), {
          status: 401,
          headers: { "Content-Type": "application/json" }
        });
      }
      
      const user = users[0];
      
      // In a real implementation, we would check the password hash
      // But for now, just return the user
      return new Response(JSON.stringify({ 
        data: { user: { id: user._id, email: user.email, role: user.role } }, 
        error: null 
      }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error signing in:', error);
      return new Response(JSON.stringify({ 
        data: null, 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  })
});

/**
 * Auth signup endpoint
 */
http.route({
  path: "/api/auth/signup",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      const body = await request.json();
      const { email, password } = body;
      
      // Check if user already exists
      const existingUsers = await ctx.db
        .query("users")
        .filter(q => q.eq(q.field("email"), email))
        .collect();
      
      if (existingUsers.length > 0) {
        return new Response(JSON.stringify({ 
          data: null, 
          error: "User already exists" 
        }), {
          status: 409,
          headers: { "Content-Type": "application/json" }
        });
      }
      
      // In a real implementation, we would hash the password
      // Create the user
      const now = Date.now();
      const userId = await ctx.db.insert("users", {
        email,
        passwordHash: password, // In reality, this would be hashed
        role: "user",
        createdAt: now,
        updatedAt: now
      });
      
      return new Response(JSON.stringify({ 
        data: { user: { id: userId, email, role: "user" } }, 
        error: null 
      }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error signing up:', error);
      return new Response(JSON.stringify({ 
        data: null, 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  })
});

/**
 * Auth signout endpoint
 */
http.route({
  path: "/api/auth/signout",
  method: "POST",
  handler: httpAction(async (ctx, request) => {
    try {
      // In a real implementation, we would invalidate the session
      return new Response(JSON.stringify({ data: {}, error: null }), {
        status: 200,
        headers: { "Content-Type": "application/json" }
      });
    } catch (error) {
      console.error('Error signing out:', error);
      return new Response(JSON.stringify({ 
        data: null, 
        error: error instanceof Error ? error.message : 'Unknown error' 
      }), {
        status: 400,
        headers: { "Content-Type": "application/json" }
      });
    }
  })
});

export default http;