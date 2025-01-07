import pytest
from fastapi.testclient import TestClient
from app.main import app  
import subprocess
import requests


# ========== Backend Tests ==========
client = TestClient(app)


def test_root_endpoint():
    """Test the root endpoint of the FastAPI backend."""
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {"message": "Welcome to DrInfo.ai"}


def test_chat_endpoint_valid():
    """Test the chat endpoint with a valid query."""
    payload = {"user_message": "treatment for hypertension"}
    response = client.post("/api/chat", json=payload)
    assert response.status_code == 200
    assert "summary" in response.json()


def test_chat_endpoint_empty_input():
    """Test the chat endpoint with empty input."""
    payload = {"user_message": ""}
    response = client.post("/api/chat", json=payload)
    assert response.status_code == 422  # Unprocessable Entity


def test_chat_endpoint_invalid_payload():
    """Test the chat endpoint with an invalid payload."""
    payload = {"invalid_key": "test"}
    response = client.post("/api/chat", json=payload)
    assert response.status_code == 422  # Unprocessable Entity


# ========== Frontend Tests ==========
@pytest.fixture
def run_streamlit():
    """Run the Streamlit app as a subprocess."""
    process = subprocess.Popen(
        ["streamlit", "run", "app/frontend.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    yield process
    process.terminate()


def test_streamlit_running(run_streamlit):
    """Check if the Streamlit app starts and is accessible."""
    import time

    time.sleep(5)  # Wait for Streamlit to start
    response = requests.get("http://localhost:8501")
    assert response.status_code == 200
    assert "DrInfo" in response.text


def test_streamlit_input(run_streamlit):
    """Check if Streamlit accepts user input."""
    response = requests.post("http://localhost:8501", json={"user_message": "test input"})
    assert response.status_code == 200


# ========== Integration Test ==========
def test_backend_and_frontend_integration(run_streamlit):
    """Test the end-to-end integration of backend and frontend."""
    payload = {"user_message": "treatment for diabetes"}
    response = requests.post("http://localhost:8501/api/chat", json=payload)
    assert response.status_code == 200
    assert "summary" in response.json()