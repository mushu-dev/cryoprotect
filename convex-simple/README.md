# Simplified Convex Implementation

This is a simplified version of the Convex implementation for CryoProtect, designed for easy deployment and testing.

## Structure

- `convex.json` - Basic Convex configuration
- `schema.js` - Simple schema definition with just a molecules table
- `hello.js` - Basic query function
- `molecules.js` - Simple CRUD functions for molecules

## Usage

Run this simplified version with:

```bash
npm run convex:simple
```

This will initialize a new Convex project and deploy this simplified implementation.

## Testing

After deployment, you can test the API with:

1. Create a molecule:
   ```javascript
   await client.mutation("molecules:createMolecule", {
     name: "Water",
     formula: "H2O"
   });
   ```

2. List all molecules:
   ```javascript
   await client.query("molecules:listMolecules");
   ```

3. Get a molecule by ID:
   ```javascript
   await client.query("molecules:getMolecule", { id: moleculeId });
   ```

## Next Steps

Once this simplified version is working correctly, you can gradually migrate the full implementation from the main `convex` directory.