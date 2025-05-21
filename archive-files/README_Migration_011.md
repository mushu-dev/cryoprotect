# Migration 011: Create predictions and experiments tables

This migration adds the missing `predictions` and `experiments` tables required for CryoProtect v2 API functionality.

## How to Apply

1. **Via Supabase Dashboard:**
   - Open the Supabase dashboard for your project.
   - Go to the SQL editor.
   - Copy and paste the contents of `migrations/011_create_predictions_experiments.sql` into the editor.
   - Run the SQL.

2. **Via psql or CLI:**
   - Connect to your database using `psql` or a compatible tool.
   - Run:
     ```
     \i migrations/011_create_predictions_experiments.sql
     ```

## Verification

After applying the migration:
- Check that the tables `public.predictions` and `public.experiments` exist.
- Run the API verification script to confirm that endpoints depending on these tables now function.

## Troubleshooting

- If you encounter permission errors, ensure you are connected as a user with schema modification rights.
- If the tables already exist, the `IF NOT EXISTS` clause will prevent errors.