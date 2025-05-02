"""
CryoProtect v2 - Backup System

This package provides a robust, production-ready backup system for CryoProtect v2.
It handles:
1. Automated backup scheduling
2. Integrity verification after each backup
3. Configurable retention policies
4. Cross-region backup storage
5. Restoration procedures with verification
"""

from .backup_manager import BackupManager

__all__ = ['BackupManager']