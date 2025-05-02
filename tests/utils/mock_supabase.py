"""
Mock Supabase Client for Testing

This module provides a mock implementation of the Supabase client
to be used in unit tests. It simulates the behavior of the real
Supabase client without making actual network requests.
"""

class MockSupabaseClient:
    """
    A mock Supabase client that simulates database operations for testing.
    """

    def __init__(self):
        """
        Initialize the mock client with in-memory storage.
        """
        self._tables = {}

    def table(self, table_name):
        """
        Return a mock table interface for the given table name.

        Args:
            table_name (str): The name of the table.

        Returns:
            MockSupabaseTable: A mock table interface.
        """
        if table_name not in self._tables:
            self._tables[table_name] = MockSupabaseTable(table_name)
        return self._tables[table_name]


class MockSupabaseTable:
    """
    A mock table interface for simulating Supabase table operations.
    """

    def __init__(self, table_name):
        """
        Initialize the mock table with an empty list of records.

        Args:
            table_name (str): The name of the table.
        """
        self.table_name = table_name
        self.records = []

    def insert(self, data):
        """
        Simulate inserting data into the table.

        Args:
            data (dict or list of dict): The data to insert.

        Returns:
            dict: A mock response indicating success.
        """
        if isinstance(data, dict):
            self.records.append(data)
        elif isinstance(data, list):
            self.records.extend(data)
        else:
            raise ValueError("Data must be a dict or list of dicts.")
        return {"data": data, "error": None}

    def select(self):
        """
        Simulate selecting all records from the table.

        Returns:
            dict: A mock response containing all records.
        """
        return {"data": self.records.copy(), "error": None}

    def delete(self, match=None):
        """
        Simulate deleting records from the table.

        Args:
            match (dict, optional): Criteria to match records for deletion.

        Returns:
            dict: A mock response indicating deleted records.
        """
        if match is None:
            deleted = self.records.copy()
            self.records.clear()
        else:
            deleted = [r for r in self.records if all(r.get(k) == v for k, v in match.items())]
            self.records = [r for r in self.records if not all(r.get(k) == v for k, v in match.items())]
        return {"data": deleted, "error": None}

    def update(self, values, match=None):
        """
        Simulate updating records in the table.

        Args:
            values (dict): The values to update.
            match (dict, optional): Criteria to match records for updating.

        Returns:
            dict: A mock response indicating updated records.
        """
        updated = []
        for record in self.records:
            if match is None or all(record.get(k) == v for k, v in match.items()):
                record.update(values)
                updated.append(record)
        return {"data": updated, "error": None}