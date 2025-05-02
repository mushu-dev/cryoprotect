#!/usr/bin/env python3
"""
Unit tests for the IP Resolution module.
"""

import unittest
from unittest.mock import patch, MagicMock
import socket
import subprocess
import platform
import os
import sys

# Add the project root directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from ip_resolver import (
    resolve_ip_address,
    _standard_dns_resolution,
    _alternative_dns_resolution,
    _resolve_with_nslookup,
    _resolve_with_dig,
    _mcp_resolution,
    _heuristic_resolution,
    update_env_with_ip,
    resolve_and_update_env
)

class TestIPResolver(unittest.TestCase):
    """Test cases for the IP Resolution module."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_hostname = "db.tsdlmynydfuypiugmkev.supabase.co"
        self.test_ip = "123.45.67.89"
        self.test_dns_server = "8.8.8.8"

    def test_standard_dns_resolution_success(self):
        """Test standard DNS resolution when successful."""
        with patch('socket.gethostbyname', return_value=self.test_ip):
            result = _standard_dns_resolution(self.test_hostname)
            self.assertEqual(result, self.test_ip)

    def test_standard_dns_resolution_failure(self):
        """Test standard DNS resolution when it fails."""
        with patch('socket.gethostbyname', side_effect=socket.gaierror):
            result = _standard_dns_resolution(self.test_hostname)
            self.assertIsNone(result)

    @patch('platform.system', return_value='Windows')
    @patch('ip_resolver._resolve_with_nslookup')
    def test_alternative_dns_resolution_windows(self, mock_nslookup, mock_platform):
        """Test alternative DNS resolution on Windows."""
        mock_nslookup.return_value = self.test_ip
        result = _alternative_dns_resolution(self.test_hostname)
        self.assertEqual(result, self.test_ip)
        mock_nslookup.assert_called_with(self.test_hostname, "8.8.8.8")

    @patch('platform.system', return_value='Linux')
    @patch('ip_resolver._resolve_with_dig')
    def test_alternative_dns_resolution_unix(self, mock_dig, mock_platform):
        """Test alternative DNS resolution on Unix."""
        mock_dig.return_value = self.test_ip
        result = _alternative_dns_resolution(self.test_hostname)
        self.assertEqual(result, self.test_ip)
        mock_dig.assert_called_with(self.test_hostname, "8.8.8.8")

    @patch('subprocess.run')
    def test_resolve_with_nslookup_success(self, mock_run):
        """Test nslookup resolution when successful."""
        # Mock the subprocess.run result
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = """
Server:  dns.google
Address:  8.8.8.8

Non-authoritative answer:
Name:    db.tsdlmynydfuypiugmkev.supabase.co
Address:  123.45.67.89
"""
        mock_run.return_value = mock_process

        result = _resolve_with_nslookup(self.test_hostname, self.test_dns_server)
        self.assertEqual(result, self.test_ip)
        mock_run.assert_called_with(
            ["nslookup", self.test_hostname, self.test_dns_server],
            capture_output=True, text=True, timeout=5
        )

    @patch('subprocess.run')
    def test_resolve_with_nslookup_failure(self, mock_run):
        """Test nslookup resolution when it fails."""
        # Mock the subprocess.run result for failure
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_run.return_value = mock_process

        result = _resolve_with_nslookup(self.test_hostname, self.test_dns_server)
        self.assertIsNone(result)

    @patch('subprocess.run')
    def test_resolve_with_dig_success(self, mock_run):
        """Test dig resolution when successful."""
        # Mock the subprocess.run result
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "123.45.67.89\n"
        mock_run.return_value = mock_process

        result = _resolve_with_dig(self.test_hostname, self.test_dns_server)
        self.assertEqual(result, self.test_ip)
        mock_run.assert_called_with(
            ["dig", "+short", self.test_hostname, f"@{self.test_dns_server}"],
            capture_output=True, text=True, timeout=5
        )

    @patch('subprocess.run')
    def test_resolve_with_dig_failure(self, mock_run):
        """Test dig resolution when it fails."""
        # Mock the subprocess.run result for failure
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_run.return_value = mock_process

        result = _resolve_with_dig(self.test_hostname, self.test_dns_server)
        self.assertIsNone(result)

    def test_mcp_resolution(self):
        """Test MCP resolution (placeholder implementation)."""
        # Since this is a placeholder, it should return None
        result = _mcp_resolution(self.test_hostname)
        self.assertIsNone(result)

    def test_heuristic_resolution(self):
        """Test heuristic resolution (placeholder implementation)."""
        # Since this is a placeholder, it should return None
        result = _heuristic_resolution(self.test_hostname)
        self.assertIsNone(result)

    @patch('os.path.exists', return_value=True)
    def test_update_env_with_ip_new_var(self, mock_exists):
        """Test updating .env file with a new IP variable."""
        # Create a mock for the open function
        mock_file = unittest.mock.mock_open(read_data="DB_HOST=localhost\n")
        
        # Apply the mock to the built-in open function
        with patch('builtins.open', mock_file):
            result = update_env_with_ip(self.test_hostname, self.test_ip)
            self.assertTrue(result)
        
        # Check that the file was written with the expected content
        # The writelines method is called with a list of lines
        write_calls = mock_file().writelines.call_args[0][0]
        self.assertIn('DB_HOST=localhost\n', write_calls)
        self.assertIn('\n# Added by IP Resolution Module\n', write_calls)
        self.assertIn(f'SUPABASE_DB_IP={self.test_ip}\n', write_calls)

    @patch('os.path.exists', return_value=True)
    def test_update_env_with_ip_existing_var(self, mock_exists):
        """Test updating .env file with an existing IP variable."""
        # Create a mock for the open function
        mock_file = unittest.mock.mock_open(read_data="SUPABASE_DB_IP=old.ip.address\n")
        
        # Apply the mock to the built-in open function
        with patch('builtins.open', mock_file):
            result = update_env_with_ip(self.test_hostname, self.test_ip)
            self.assertTrue(result)
        
        # Check that the file was written with the expected content
        # The writelines method is called with a list of lines
        write_calls = mock_file().writelines.call_args[0][0]
        self.assertIn(f'SUPABASE_DB_IP={self.test_ip}\n', write_calls)
        self.assertNotIn('old.ip.address', str(write_calls))

    @patch('ip_resolver.resolve_ip_address', return_value=None)
    def test_resolve_and_update_env_failure(self, mock_resolve):
        """Test resolve_and_update_env when resolution fails."""
        result = resolve_and_update_env(self.test_hostname)
        self.assertEqual(result['hostname'], self.test_hostname)
        self.assertIsNone(result['ip'])
        self.assertFalse(result['updated'])

    @patch('ip_resolver.resolve_ip_address', return_value="123.45.67.89")
    @patch('ip_resolver.update_env_with_ip', return_value=True)
    def test_resolve_and_update_env_success(self, mock_update, mock_resolve):
        """Test resolve_and_update_env when successful."""
        result = resolve_and_update_env(self.test_hostname)
        self.assertEqual(result['hostname'], self.test_hostname)
        self.assertEqual(result['ip'], self.test_ip)
        self.assertTrue(result['updated'])
        mock_update.assert_called_with(self.test_hostname, self.test_ip)

    @patch('ip_resolver._standard_dns_resolution', return_value=None)
    @patch('ip_resolver._alternative_dns_resolution', return_value=None)
    @patch('ip_resolver._mcp_resolution', return_value=None)
    @patch('ip_resolver._heuristic_resolution', return_value=None)
    def test_resolve_ip_address_all_methods_fail(self, mock_heuristic, mock_mcp, mock_alt, mock_std):
        """Test resolve_ip_address when all methods fail."""
        result = resolve_ip_address(self.test_hostname)
        self.assertIsNone(result)
        mock_std.assert_called_once_with(self.test_hostname)
        mock_alt.assert_called_once_with(self.test_hostname)
        mock_mcp.assert_called_once_with(self.test_hostname)
        mock_heuristic.assert_called_once_with(self.test_hostname)

    @patch('ip_resolver._standard_dns_resolution', return_value="123.45.67.89")
    @patch('ip_resolver._alternative_dns_resolution')
    def test_resolve_ip_address_first_method_succeeds(self, mock_alt, mock_std):
        """Test resolve_ip_address when the first method succeeds."""
        result = resolve_ip_address(self.test_hostname)
        self.assertEqual(result, self.test_ip)
        mock_std.assert_called_once_with(self.test_hostname)
        mock_alt.assert_not_called()

if __name__ == '__main__':
    unittest.main()