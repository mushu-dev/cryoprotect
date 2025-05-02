"""
test_codex_cli_wrapper.py

Automated tests for codex_cli_wrapper.py, covering integration with the OpenAI Codex CLI via subprocess.
All subprocess and environment variable interactions are mocked to ensure tests are self-contained and do not require the actual Codex CLI or API key.

Tested Scenarios:
1. Successful invocation of run_codex (mock subprocess to return expected output).
2. Error when OPENAI_API_KEY is missing (should raise EnvironmentError).
3. Error when Codex CLI is not installed (mock FileNotFoundError, should raise CodexCLIError).
4. Error when Codex CLI returns nonzero exit code (mock subprocess to return error, should raise CodexCLIError).
"""

import os
import pytest
from unittest import mock

from codex_cli_wrapper import run_codex, CodexCLIError

@pytest.fixture
def mock_env_with_api_key(monkeypatch):
    """Fixture to set OPENAI_API_KEY in the environment."""
    monkeypatch.setenv("OPENAI_API_KEY", "test-api-key")
    yield
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)

def test_run_codex_success(mock_env_with_api_key):
    """
    Test successful invocation of run_codex.
    Mocks subprocess.run to return expected output.
    """
    expected_output = "Codex result here"
    mock_result = mock.Mock()
    mock_result.returncode = 0
    mock_result.stdout = expected_output
    mock_result.stderr = ""

    with mock.patch("subprocess.run", return_value=mock_result) as mock_run:
        result = run_codex("Test prompt", output_format="text")
        assert result == expected_output
        # Check CLI command construction
        mock_run.assert_called_once()
        args, kwargs = mock_run.call_args
        assert "codex" in args[0]
        assert "Test prompt" in args[0]

def test_run_codex_success_json(mock_env_with_api_key):
    """
    Test successful invocation of run_codex with output_format='json'.
    Ensures '--json' is passed to the CLI.
    """
    expected_output = '{"result": "Codex JSON"}'
    mock_result = mock.Mock()
    mock_result.returncode = 0
    mock_result.stdout = expected_output
    mock_result.stderr = ""

    with mock.patch("subprocess.run", return_value=mock_result) as mock_run:
        result = run_codex("Test prompt", output_format="json")
        assert result == expected_output
        args, kwargs = mock_run.call_args
        assert "--json" in args[0]

def test_run_codex_missing_api_key(monkeypatch):
    """
    Test error when OPENAI_API_KEY is missing.
    Should raise EnvironmentError.
    """
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    with pytest.raises(EnvironmentError) as excinfo:
        run_codex("Test prompt")
    assert "OPENAI_API_KEY environment variable is not set" in str(excinfo.value)

def test_run_codex_cli_not_installed(mock_env_with_api_key):
    """
    Test error when Codex CLI is not installed.
    Mocks subprocess.run to raise FileNotFoundError.
    Should raise CodexCLIError.
    """
    with mock.patch("subprocess.run", side_effect=FileNotFoundError("No such file or directory: 'codex'")):
        with pytest.raises(CodexCLIError) as excinfo:
            run_codex("Test prompt")
    assert "Codex CLI is not installed or not found in your PATH" in str(excinfo.value)

def test_run_codex_cli_failure(mock_env_with_api_key):
    """
    Test error when Codex CLI returns nonzero exit code.
    Mocks subprocess.run to return error.
    Should raise CodexCLIError with stderr and stdout in the message.
    """
    mock_result = mock.Mock()
    mock_result.returncode = 1
    mock_result.stdout = "Some output"
    mock_result.stderr = "Error: something went wrong"

    with mock.patch("subprocess.run", return_value=mock_result):
        with pytest.raises(CodexCLIError) as excinfo:
            run_codex("Test prompt")
    msg = str(excinfo.value)
    assert "Codex CLI failed with exit code 1" in msg
    assert "Error: something went wrong" in msg
    assert "Some output" in msg