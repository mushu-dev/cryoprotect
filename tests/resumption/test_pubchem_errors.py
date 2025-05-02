import pytest

class MockPubChemImporter:
    def fetch_properties(self, cid):
        # Simulate different error scenarios based on CID
        if cid == "503":
            raise ConnectionError("PubChem API returned 503 Service Unavailable")
        elif cid == "bad_data":
            raise ValueError("Error in molecule data")
        elif cid == "out_of_range":
            return {"LogP": 7.0, "TPSA": 100, "MW": 500}
        elif cid == "none_value":
            return {"LogP": None, "TPSA": 100, "MW": 500}
        else:
            return {"LogP": 2.0, "TPSA": 100, "MW": 500}

def validate_properties(props):
    # Simulate property validation logic
    if props["LogP"] is None:
        raise ValueError("LogP None outside range (-5, 5)")
    if not (-5 <= props["LogP"] <= 5):
        raise ValueError(f"LogP {props['LogP']} outside range (-5, 5)")
    return True

def test_api_503_error():
    """
    Test handling of PubChem API 503 Service Unavailable error.
    """
    importer = MockPubChemImporter()
    with pytest.raises(ConnectionError, match="503 Service Unavailable"):
        importer.fetch_properties("503")

def test_molecule_data_error():
    """
    Test handling of molecule data parsing or missing/invalid fields.
    """
    importer = MockPubChemImporter()
    with pytest.raises(ValueError, match="Error in molecule data"):
        importer.fetch_properties("bad_data")

def test_property_out_of_range_error():
    """
    Test handling of property out-of-range error (e.g., LogP 7.0 outside range).
    """
    importer = MockPubChemImporter()
    props = importer.fetch_properties("out_of_range")
    with pytest.raises(ValueError, match="outside range"):
        validate_properties(props)

def test_none_value_error():
    """
    Test handling of None value for a property (e.g., LogP None outside range).
    """
    importer = MockPubChemImporter()
    props = importer.fetch_properties("none_value")
    with pytest.raises(ValueError, match="None outside range"):
        validate_properties(props)