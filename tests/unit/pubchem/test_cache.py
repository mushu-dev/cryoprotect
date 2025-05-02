import os
import sqlite3
import tempfile
import unittest
from pathlib import Path

from pubchem.cache import initialize_cache, get_compound, get_cache_stats
import json
import time


class TestPubChemCacheInitialization(unittest.TestCase):
    """Test cases for PubChem cache initialization and schema creation."""

    def setUp(self):
        """Create a temporary directory for test database files."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.temp_dir.name, "test_pubchem_cache.sqlite")

    def tearDown(self):
        """Clean up temporary files."""
        self.temp_dir.cleanup()

    def test_initialize_cache_creates_database(self):
        """Test that initialize_cache creates the database file if it doesn't exist."""
        # Ensure the file doesn't exist
        if os.path.exists(self.db_path):
            os.remove(self.db_path)

        # Initialize the cache
        initialize_cache(self.db_path)

        # Check that the file was created
        self.assertTrue(os.path.exists(self.db_path))
        self.assertTrue(os.path.isfile(self.db_path))

    def test_initialize_cache_creates_tables(self):
        """Test that initialize_cache creates the required tables."""
        # Initialize the cache
        initialize_cache(self.db_path)

        # Connect to the database and check tables
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get list of tables
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]

        # Check that required tables exist
        self.assertIn("compound_cache", tables)
        self.assertIn("cache_stats", tables)

        # Close connection
        conn.close()

    def test_compound_cache_schema(self):
        """Test that the compound_cache table has the correct schema."""
        # Initialize the cache
        initialize_cache(self.db_path)

        # Connect to the database and check schema
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get schema for compound_cache table
        cursor.execute("PRAGMA table_info(compound_cache)")
        columns = {row[1]: row[2] for row in cursor.fetchall()}

        # Check columns and types
        self.assertIn("cid", columns)
        self.assertIn("data", columns)
        self.assertIn("last_updated", columns)
        self.assertIn("last_accessed", columns)
        self.assertIn("source", columns)

        # Check column types
        self.assertEqual(columns["cid"], "INTEGER")
        self.assertEqual(columns["data"], "TEXT")
        self.assertEqual(columns["last_updated"], "INTEGER")
        self.assertEqual(columns["last_accessed"], "INTEGER")
        self.assertEqual(columns["source"], "TEXT")

        # Close connection
        conn.close()

    def test_cache_stats_schema(self):
        """Test that the cache_stats table has the correct schema."""
        # Initialize the cache
        initialize_cache(self.db_path)

        # Connect to the database and check schema
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get schema for cache_stats table
        cursor.execute("PRAGMA table_info(cache_stats)")
        columns = {row[1]: row[2] for row in cursor.fetchall()}

        # Check columns and types
        self.assertIn("key", columns)
        self.assertIn("value", columns)

        # Check column types
        self.assertEqual(columns["key"], "TEXT")
        self.assertEqual(columns["value"], "INTEGER")

        # Close connection
        conn.close()

    def test_indexes_created(self):
        """Test that the required indexes are created."""
        # Initialize the cache
        initialize_cache(self.db_path)

        # Connect to the database and check indexes
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get indexes
        cursor.execute("SELECT name FROM sqlite_master WHERE type='index'")
        indexes = [row[0] for row in cursor.fetchall()]

        # Check that required indexes exist
        self.assertIn("idx_last_accessed", indexes)

        # Close connection
        conn.close()

    def test_stats_initialized(self):
        """Test that the cache stats are initialized correctly."""
        # Initialize the cache
        initialize_cache(self.db_path)

        # Connect to the database and check stats
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get stats
        cursor.execute("SELECT key, value FROM cache_stats")
        stats = {row[0]: row[1] for row in cursor.fetchall()}

        # Check that required stats exist and are initialized to 0
        self.assertIn("hits", stats)
        self.assertIn("misses", stats)
        self.assertIn("total_entries", stats)
        self.assertIn("hit_rate", stats)
        self.assertEqual(stats["hits"], 0)
        self.assertEqual(stats["misses"], 0)
        self.assertEqual(stats["total_entries"], 0)
        self.assertEqual(stats["hit_rate"], 0)

        # Close connection
        conn.close()

    def test_repeated_initialization(self):
        """Test that repeated initialization doesn't cause errors."""
        # Initialize the cache multiple times
        initialize_cache(self.db_path)
        initialize_cache(self.db_path)
        initialize_cache(self.db_path)

        # Connect to the database and check that it's still valid
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Check that we can query the database
        cursor.execute("SELECT count(*) FROM sqlite_master")
        count = cursor.fetchone()[0]
        self.assertGreater(count, 0)

        # Close connection
        conn.close()

    def test_default_path(self):
        """Test initialization with default path."""
        # Use a temporary directory as the current directory
        original_dir = os.getcwd()
        os.chdir(self.temp_dir.name)

        try:
            # Initialize with default path
            initialize_cache()

            # Check that the file was created in the expected location
            expected_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                        '..', '..', '..', 'pubchem', 'pubchem_cache.sqlite')
            self.assertTrue(os.path.exists(expected_path) or 
                           os.path.exists(os.path.join(self.temp_dir.name, 'pubchem_cache.sqlite')))
        finally:
            # Restore original directory
            os.chdir(original_dir)


class TestPubChemCacheRetrieval(unittest.TestCase):
    """Test cases for PubChem cache retrieval functionality."""

    def setUp(self):
        """Create a temporary directory and initialize the cache with test data."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.temp_dir.name, "test_pubchem_cache.sqlite")
        
        # Initialize the cache
        initialize_cache(self.db_path)
        
        # Insert test data directly into the database
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
        # Current time for testing
        self.current_time = int(time.time())
        
        # Insert test compounds
        test_compounds = [
            (1, json.dumps({"name": "Compound 1", "formula": "C1"}), self.current_time, self.current_time, "pubchem"),
            (2, json.dumps({"name": "Compound 2", "formula": "C2"}), self.current_time - 3600, self.current_time - 3600, "pubchem"),
            (3, json.dumps({"name": "Compound 3", "formula": "C3"}), self.current_time - 86400, self.current_time - 86400, "pubchem"),
        ]
        
        for compound in test_compounds:
            self.cursor.execute('''
            INSERT INTO compound_cache (cid, data, last_updated, last_accessed, source)
            VALUES (?, ?, ?, ?, ?)
            ''', compound)
        
        # Reset stats
        self.cursor.execute("UPDATE cache_stats SET value = 0 WHERE key IN ('hits', 'misses', 'hit_rate')")
        self.cursor.execute("UPDATE cache_stats SET value = 3 WHERE key = 'total_entries'")
        
        self.conn.commit()
        
        # Monkey patch the module to use our test database
        self._original_dirname = os.path.dirname
        os.path.dirname = lambda x: self.temp_dir.name if x == __file__ else self._original_dirname(x)

    def tearDown(self):
        """Clean up temporary files and restore original functions."""
        self.conn.close()
        self.temp_dir.cleanup()
        os.path.dirname = self._original_dirname

    def test_get_existing_compound(self):
        """Test retrieving a compound that exists in the cache."""
        # Get compound with ID 1
        result = get_compound(1, db_path=self.db_path)
        
        # Check that the result is correct
        self.assertIsNotNone(result)
        self.assertEqual(result["name"], "Compound 1")
        self.assertEqual(result["formula"], "C1")
        
        # Get cache stats and check values
        stats = get_cache_stats(self.db_path)
        
        # Check that hit count was incremented
        self.assertEqual(stats['hits'], 1)
        
        # Check that miss count was not incremented
        self.assertEqual(stats['misses'], 0)
        
        # Check that hit rate was updated
        self.assertEqual(stats['hit_rate'], 100)  # 1 hit, 0 misses = 100% hit rate

    def test_get_nonexistent_compound(self):
        """Test retrieving a compound that doesn't exist in the cache."""
        # Get compound with ID 999 (doesn't exist)
        result = get_compound(999, db_path=self.db_path)
        
        # Check that the result is None
        self.assertIsNone(result)
        
        # Get cache stats and check values
        stats = get_cache_stats(self.db_path)
        
        # Check that hit count was not incremented
        self.assertEqual(stats['hits'], 0)
        
        # Check that miss count was incremented
        self.assertEqual(stats['misses'], 1)
        
        # Check that hit rate was updated
        self.assertEqual(stats['hit_rate'], 0)  # 0 hits, 1 miss = 0% hit rate

    def test_get_expired_compound(self):
        """Test retrieving a compound that exists but is expired."""
        # Get compound with ID 3 with max_age_seconds of 3600 (1 hour)
        # The compound was last updated 1 day ago, so it should be considered expired
        result = get_compound(3, max_age_seconds=3600, db_path=self.db_path)
        
        # Check that the result is None (expired)
        self.assertIsNone(result)
        
        # Get cache stats and check values
        stats = get_cache_stats(self.db_path)
        
        # Check that hit count was not incremented
        self.assertEqual(stats['hits'], 0)
        
        # Check that miss count was incremented
        self.assertEqual(stats['misses'], 1)

    def test_get_recent_compound_with_expiry(self):
        """Test retrieving a compound that exists and is not expired."""
        # Get compound with ID 1 with max_age_seconds of 3600 (1 hour)
        # The compound was just updated, so it should not be expired
        result = get_compound(1, max_age_seconds=3600, db_path=self.db_path)
        
        # Check that the result is correct
        self.assertIsNotNone(result)
        self.assertEqual(result["name"], "Compound 1")
        
        # Get cache stats and check values
        stats = get_cache_stats(self.db_path)
        
        # Check that hit count was incremented
        self.assertEqual(stats['hits'], 1)

    def test_last_accessed_updated(self):
        """Test that the last_accessed timestamp is updated on retrieval."""
        # Get the initial last_accessed timestamp
        self.cursor.execute("SELECT last_accessed FROM compound_cache WHERE cid = 2")
        initial_last_accessed = self.cursor.fetchone()[0]
        
        # Wait a short time to ensure timestamp will be different
        time.sleep(0.1)
        
        # Get compound with ID 2
        result = get_compound(2, db_path=self.db_path)
        
        # Check that the result is correct
        self.assertIsNotNone(result)
        
        # Get the updated last_accessed timestamp
        self.cursor.execute("SELECT last_accessed FROM compound_cache WHERE cid = 2")
        updated_last_accessed = self.cursor.fetchone()[0]
        
        # Check that the timestamp was updated
        self.assertGreater(updated_last_accessed, initial_last_accessed)

    def test_multiple_hits_and_misses(self):
        """Test that hit/miss statistics are updated correctly with multiple operations."""
        # Get existing compounds
        get_compound(1, db_path=self.db_path)  # Hit
        get_compound(2, db_path=self.db_path)  # Hit
        
        # Get non-existent compounds
        get_compound(999, db_path=self.db_path)  # Miss
        get_compound(998, db_path=self.db_path)  # Miss
        
        # Get expired compound
        get_compound(3, max_age_seconds=3600, db_path=self.db_path)  # Miss (expired)
        
        # Get cache stats and check values
        stats = get_cache_stats(self.db_path)
        
        # Check hit count
        self.assertEqual(stats['hits'], 2)
        
        # Check miss count
        self.assertEqual(stats['misses'], 3)
        
        # Check hit rate
        self.assertEqual(stats['hit_rate'], 40)  # 2 hits, 3 misses = 40% hit rate


class TestPubChemCacheStorage(unittest.TestCase):
    """Test cases for PubChem cache storage functionality."""

    def setUp(self):
        """Create a temporary directory and initialize the cache."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.temp_dir.name, "test_pubchem_cache.sqlite")
        
        # Initialize the cache
        initialize_cache(self.db_path)
        
        # Connect to the database for direct verification
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
        # Current time for testing
        self.current_time = int(time.time())

    def tearDown(self):
        """Clean up temporary files."""
        self.conn.close()
        self.temp_dir.cleanup()

    def test_store_new_compound(self):
        """Test storing a new compound that doesn't exist in the cache."""
        # Test data
        cid = 12345
        data = {"name": "Test Compound", "formula": "C10H15N5O10P2", "molecular_weight": 507.18}
        
        # Store the compound
        from pubchem.cache import store_compound
        store_compound(cid, data, db_path=self.db_path)
        
        # Verify the compound was stored
        self.cursor.execute("SELECT cid, data, source FROM compound_cache WHERE cid = ?", (cid,))
        result = self.cursor.fetchone()
        
        # Check that the compound exists
        self.assertIsNotNone(result)
        self.assertEqual(result[0], cid)
        
        # Check that the data was stored correctly
        stored_data = json.loads(result[1])
        self.assertEqual(stored_data, data)
        
        # Check that the source is correct
        self.assertEqual(result[2], "pubchem")
        
        # Get cache stats and check total entries
        from pubchem.cache import get_cache_stats
        stats = get_cache_stats(self.db_path)
        self.assertEqual(stats["total_entries"], 1)

    def test_update_existing_compound(self):
        """Test updating a compound that already exists in the cache."""
        # Test data
        cid = 54321
        initial_data = {"name": "Initial Compound", "formula": "C8H10N4O2"}
        updated_data = {"name": "Updated Compound", "formula": "C8H10N4O2", "molecular_weight": 194.19}
        
        # Store the initial compound
        from pubchem.cache import store_compound
        store_compound(cid, initial_data, db_path=self.db_path)
        
        # Get the initial timestamp
        self.cursor.execute("SELECT last_updated FROM compound_cache WHERE cid = ?", (cid,))
        initial_timestamp = self.cursor.fetchone()[0]
        
        # Wait a sufficient time to ensure timestamp will be different
        # Since timestamps are stored as integers (seconds), we need to wait at least 1 second
        time.sleep(1.1)
        
        # Update the compound
        store_compound(cid, updated_data, db_path=self.db_path)
        
        # Verify the compound was updated
        self.cursor.execute("SELECT data, last_updated FROM compound_cache WHERE cid = ?", (cid,))
        result = self.cursor.fetchone()
        
        # Check that the data was updated
        stored_data = json.loads(result[0])
        self.assertEqual(stored_data, updated_data)
        
        # Check that the timestamp was updated
        self.assertGreater(result[1], initial_timestamp)
        
        # Get cache stats and check total entries (should still be 1)
        from pubchem.cache import get_cache_stats
        stats = get_cache_stats(self.db_path)
        self.assertEqual(stats["total_entries"], 1)

    def test_store_with_custom_source(self):
        """Test storing a compound with a custom source."""
        # Test data
        cid = 98765
        data = {"name": "RDKit Compound", "formula": "C6H6"}
        source = "rdkit"
        
        # Store the compound with custom source
        from pubchem.cache import store_compound
        store_compound(cid, data, source=source, db_path=self.db_path)
        
        # Verify the source was stored correctly
        self.cursor.execute("SELECT source FROM compound_cache WHERE cid = ?", (cid,))
        result = self.cursor.fetchone()
        
        # Check that the source is correct
        self.assertEqual(result[0], source)

    def test_store_multiple_compounds(self):
        """Test storing multiple compounds."""
        # Test data
        compounds = [
            (1001, {"name": "Compound 1", "formula": "C1"}),
            (1002, {"name": "Compound 2", "formula": "C2"}),
            (1003, {"name": "Compound 3", "formula": "C3"})
        ]
        
        # Store all compounds
        from pubchem.cache import store_compound
        for cid, data in compounds:
            store_compound(cid, data, db_path=self.db_path)
        
        # Verify all compounds were stored
        self.cursor.execute("SELECT COUNT(*) FROM compound_cache")
        count = self.cursor.fetchone()[0]
        self.assertEqual(count, len(compounds))
        
        # Get cache stats and check total entries
        from pubchem.cache import get_cache_stats
        stats = get_cache_stats(self.db_path)
        self.assertEqual(stats["total_entries"], len(compounds))


class TestPubChemCacheStatistics(unittest.TestCase):
    """Test cases for PubChem cache statistics reporting functionality."""

    def setUp(self):
        """Create a temporary directory and initialize the cache with test data."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.temp_dir.name, "test_pubchem_cache.sqlite")
        
        # Initialize the cache
        initialize_cache(self.db_path)
        
        # Insert test data directly into the database
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
        # Current time for testing
        self.current_time = int(time.time())
        
        # Insert test compounds with different sources and timestamps
        test_compounds = [
            (1, json.dumps({"name": "Compound 1", "formula": "C1"}), self.current_time, self.current_time, "pubchem"),
            (2, json.dumps({"name": "Compound 2", "formula": "C2"}), self.current_time - 3600, self.current_time - 1800, "pubchem"),
            (3, json.dumps({"name": "Compound 3", "formula": "C3"}), self.current_time - 86400, self.current_time - 43200, "pubchem"),
            (4, json.dumps({"name": "Compound 4", "formula": "C4"}), self.current_time - 7200, self.current_time - 3600, "rdkit"),
            (5, json.dumps({"name": "Compound 5", "formula": "C5"}), self.current_time - 14400, self.current_time - 7200, "rdkit"),
        ]
        
        for compound in test_compounds:
            self.cursor.execute('''
            INSERT INTO compound_cache (cid, data, last_updated, last_accessed, source)
            VALUES (?, ?, ?, ?, ?)
            ''', compound)
        
        # Set some test statistics
        self.cursor.execute("UPDATE cache_stats SET value = 10 WHERE key = 'hits'")
        self.cursor.execute("UPDATE cache_stats SET value = 5 WHERE key = 'misses'")
        self.cursor.execute("UPDATE cache_stats SET value = 5 WHERE key = 'total_entries'")
        self.cursor.execute("UPDATE cache_stats SET value = 67 WHERE key = 'hit_rate'")
        
        self.conn.commit()

    def tearDown(self):
        """Clean up temporary files."""
        self.conn.close()
        self.temp_dir.cleanup()

    def test_get_cache_stats(self):
        """Test that get_cache_stats returns the correct basic statistics."""
        from pubchem.cache import get_cache_stats
        
        # Get cache stats
        stats = get_cache_stats(self.db_path)
        
        # Check basic stats
        self.assertEqual(stats['hits'], 10)
        self.assertEqual(stats['misses'], 5)
        self.assertEqual(stats['total_entries'], 5)  # Should be updated to match actual count
        self.assertEqual(stats['hit_rate'], 66)  # Recalculated: int(10/(10+5)*100) = 66
    
    def test_get_cache_stats_detailed(self):
        """Test that get_cache_stats_detailed returns the correct detailed statistics."""
        from pubchem.cache import get_cache_stats_detailed
        
        # Get detailed cache stats
        detailed_stats = get_cache_stats_detailed(self.db_path)
        
        # Check basic stats are included
        self.assertEqual(detailed_stats['hits'], 10)
        self.assertEqual(detailed_stats['misses'], 5)
        self.assertEqual(detailed_stats['total_entries'], 5)
        self.assertEqual(detailed_stats['hit_rate'], 66)
        
        # Check age statistics
        self.assertIn('avg_age_seconds', detailed_stats)
        self.assertIn('min_age_seconds', detailed_stats)
        self.assertIn('max_age_seconds', detailed_stats)
        self.assertIn('avg_access_age_seconds', detailed_stats)
        self.assertIn('min_access_age_seconds', detailed_stats)
        self.assertIn('max_access_age_seconds', detailed_stats)
        
        # Verify min/max age values
        self.assertEqual(detailed_stats['min_age_seconds'], 0)  # Compound 1 was just updated
        self.assertGreaterEqual(detailed_stats['max_age_seconds'], 86400)  # Compound 3 is at least 1 day old
        
        # Check source distribution
        self.assertIn('source_distribution', detailed_stats)
        self.assertEqual(detailed_stats['source_distribution']['pubchem'], 3)
        self.assertEqual(detailed_stats['source_distribution']['rdkit'], 2)
        
        # Check size statistics
        self.assertIn('total_data_size_bytes', detailed_stats)
        self.assertIn('avg_data_size_bytes', detailed_stats)
        self.assertIn('min_data_size_bytes', detailed_stats)
        self.assertIn('max_data_size_bytes', detailed_stats)
        self.assertGreater(detailed_stats['total_data_size_bytes'], 0)
        
        # Check recent/oldest access lists
        self.assertIn('recent_access', detailed_stats)
        self.assertIn('oldest_access', detailed_stats)
        self.assertEqual(len(detailed_stats['recent_access']), 5)
        self.assertEqual(len(detailed_stats['oldest_access']), 5)
        
        # Most recently accessed should be Compound 1
        self.assertEqual(detailed_stats['recent_access'][0]['cid'], 1)
        
        # Least recently accessed should be Compound 3
        self.assertEqual(detailed_stats['oldest_access'][0]['cid'], 3)
        
        # Check database file size
        self.assertIn('db_file_size_bytes', detailed_stats)
        self.assertGreater(detailed_stats['db_file_size_bytes'], 0)

    def test_empty_cache_stats(self):
        """Test that get_cache_stats_detailed handles an empty cache gracefully."""
        # Create a new empty cache
        empty_db_path = os.path.join(self.temp_dir.name, "empty_cache.sqlite")
        initialize_cache(empty_db_path)
        
        from pubchem.cache import get_cache_stats_detailed
        
        # Get detailed cache stats for empty cache
        detailed_stats = get_cache_stats_detailed(empty_db_path)
        
        # Check basic stats
        self.assertEqual(detailed_stats['hits'], 0)
        self.assertEqual(detailed_stats['misses'], 0)
        self.assertEqual(detailed_stats['total_entries'], 0)
        self.assertEqual(detailed_stats['hit_rate'], 0)
        
        # Check source distribution is empty
        self.assertEqual(detailed_stats['source_distribution'], {})
        
        # Check recent/oldest access lists are empty
        self.assertEqual(detailed_stats['recent_access'], [])
        self.assertEqual(detailed_stats['oldest_access'], [])


class TestPubChemCachePruning(unittest.TestCase):
    """Test cases for PubChem cache pruning functionality."""

    def setUp(self):
        """Create a temporary directory and initialize the cache with test data."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.temp_dir.name, "test_pubchem_cache.sqlite")
        
        # Initialize the cache
        initialize_cache(self.db_path)
        
        # Insert test data directly into the database
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
        # Current time for testing
        self.current_time = int(time.time())
        
        # Insert test compounds with different access timestamps
        # Compounds are ordered by last_accessed from oldest to newest
        test_compounds = [
            (1, json.dumps({"name": "Compound 1"}), self.current_time - 86400, self.current_time - 86400, "pubchem"),  # Oldest
            (2, json.dumps({"name": "Compound 2"}), self.current_time - 43200, self.current_time - 43200, "pubchem"),
            (3, json.dumps({"name": "Compound 3"}), self.current_time - 21600, self.current_time - 21600, "pubchem"),
            (4, json.dumps({"name": "Compound 4"}), self.current_time - 10800, self.current_time - 10800, "pubchem"),
            (5, json.dumps({"name": "Compound 5"}), self.current_time - 3600, self.current_time - 3600, "pubchem"),
            (6, json.dumps({"name": "Compound 6"}), self.current_time - 1800, self.current_time - 1800, "pubchem"),
            (7, json.dumps({"name": "Compound 7"}), self.current_time - 900, self.current_time - 900, "pubchem"),
            (8, json.dumps({"name": "Compound 8"}), self.current_time - 300, self.current_time - 300, "pubchem"),
            (9, json.dumps({"name": "Compound 9"}), self.current_time - 60, self.current_time - 60, "pubchem"),
            (10, json.dumps({"name": "Compound 10"}), self.current_time, self.current_time, "pubchem"),  # Newest
        ]
        
        for compound in test_compounds:
            self.cursor.execute('''
            INSERT INTO compound_cache (cid, data, last_updated, last_accessed, source)
            VALUES (?, ?, ?, ?, ?)
            ''', compound)
        
        self.conn.commit()

    def tearDown(self):
        """Clean up temporary files."""
        self.conn.close()
        self.temp_dir.cleanup()

    def test_prune_cache_no_pruning_needed(self):
        """Test that prune_cache does nothing when the cache size is below the limit."""
        from pubchem.cache import prune_cache
        
        # Prune cache with max_entries = 15 (more than the 10 we have)
        removed = prune_cache(15, self.db_path)
        
        # Check that no entries were removed
        self.assertEqual(removed, 0)
        
        # Verify that all compounds are still in the database
        self.cursor.execute("SELECT COUNT(*) FROM compound_cache")
        count = self.cursor.fetchone()[0]
        self.assertEqual(count, 10)

    def test_prune_cache_removes_oldest(self):
        """Test that prune_cache removes the least recently accessed entries."""
        from pubchem.cache import prune_cache
        
        # Prune cache with max_entries = 7 (should remove 3 oldest entries)
        removed = prune_cache(7, self.db_path)
        
        # Check that 3 entries were removed
        self.assertEqual(removed, 3)
        
        # Verify that the total count is now 7
        self.cursor.execute("SELECT COUNT(*) FROM compound_cache")
        count = self.cursor.fetchone()[0]
        self.assertEqual(count, 7)
        
        # Verify that the oldest compounds (1, 2, 3) were removed
        for cid in [1, 2, 3]:
            self.cursor.execute("SELECT COUNT(*) FROM compound_cache WHERE cid = ?", (cid,))
            count = self.cursor.fetchone()[0]
            self.assertEqual(count, 0, f"Compound {cid} should have been removed")
        
        # Verify that the newer compounds (4-10) are still present
        for cid in range(4, 11):
            self.cursor.execute("SELECT COUNT(*) FROM compound_cache WHERE cid = ?", (cid,))
            count = self.cursor.fetchone()[0]
            self.assertEqual(count, 1, f"Compound {cid} should still be present")

    def test_prune_cache_removes_all_but_max(self):
        """Test that prune_cache correctly handles the case where max_entries is very small."""
        from pubchem.cache import prune_cache
        
        # Prune cache with max_entries = 2 (should keep only the 2 newest entries)
        removed = prune_cache(2, self.db_path)
        
        # Check that 8 entries were removed
        self.assertEqual(removed, 8)
        
        # Verify that the total count is now 2
        self.cursor.execute("SELECT COUNT(*) FROM compound_cache")
        count = self.cursor.fetchone()[0]
        self.assertEqual(count, 2)
        
        # Verify that only the newest compounds (9, 10) remain
        for cid in range(1, 9):
            self.cursor.execute("SELECT COUNT(*) FROM compound_cache WHERE cid = ?", (cid,))
            count = self.cursor.fetchone()[0]
            self.assertEqual(count, 0, f"Compound {cid} should have been removed")
        
        for cid in [9, 10]:
            self.cursor.execute("SELECT COUNT(*) FROM compound_cache WHERE cid = ?", (cid,))
            count = self.cursor.fetchone()[0]
            self.assertEqual(count, 1, f"Compound {cid} should still be present")

    def test_prune_cache_empty_cache(self):
        """Test that prune_cache handles an empty cache gracefully."""
        # Create a new empty cache
        empty_db_path = os.path.join(self.temp_dir.name, "empty_cache.sqlite")
        initialize_cache(empty_db_path)
        
        from pubchem.cache import prune_cache
        
        # Prune the empty cache
        removed = prune_cache(5, empty_db_path)
        
        # Check that no entries were removed
        self.assertEqual(removed, 0)


if __name__ == "__main__":
    unittest.main()