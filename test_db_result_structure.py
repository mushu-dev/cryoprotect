from database.utils.connection import create_connection

conn = create_connection()
result = conn.sql("SELECT tablename FROM pg_tables WHERE schemaname='public' ORDER BY tablename").execute()
print("type(result.data):", type(result.data))
print("result.data:", result.data)