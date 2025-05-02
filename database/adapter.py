from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Union, Tuple

class DatabaseAdapter(ABC):
    """Abstract database adapter interface."""
    
    @abstractmethod
    def connect(self) -> bool:
        """Establish connection to the database."""
        pass
        
    @abstractmethod
    def disconnect(self) -> bool:
        """Close database connection."""
        pass
        
    @abstractmethod
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """Execute SQL query and return results."""
        pass
        
    @abstractmethod
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """Execute multiple SQL queries and return results."""
        pass
        
    @abstractmethod
    def begin_transaction(self) -> Any:
        """Begin a database transaction."""
        pass
        
    @abstractmethod
    def commit_transaction(self, transaction: Any) -> bool:
        """Commit a database transaction."""
        pass
        
    @abstractmethod
    def rollback_transaction(self, transaction: Any) -> bool:
        """Rollback a database transaction."""
        pass
        
    @abstractmethod
    def get_connection_info(self) -> Dict[str, Any]:
        """Get connection information."""
        pass
        
    @abstractmethod
    def test_connection(self) -> Tuple[bool, str]:
        """Test database connection and return status with message."""
        pass