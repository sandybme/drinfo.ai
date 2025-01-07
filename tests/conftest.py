import pytest


@pytest.fixture(scope="session")
def sample_payload():
    """Provide a sample payload for testing."""
    return {"user_message": "treatment for hypertension"}