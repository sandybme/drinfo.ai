"""
streamlit_app.py

This module contains the Streamlit-based frontend for the DrInfo.ai application. It allows users
to submit clinical queries and displays summarized responses retrieved from a backend API.

Modules:
- `requests`: Handles API requests to the backend.
- `dotenv`: Loads environment variables for configuration.
- `streamlit`: Provides the frontend framework.

Environment Variables:
- API_URL: The base URL of the FastAPI backend.

Author:
- Your Sandhanakrishnan Ravichandran
- Date: 2025-01-05
"""

import streamlit as st
import requests
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Backend API URL (set in .env)
API_URL = os.getenv("API_URL")

# Streamlit App Configuration
st.set_page_config(
    page_title="DrInfo.ai",
    layout="wide",
    page_icon="🩺",
)

# Initialize session state for chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# App Title and Description
st.title("DrInfo.ai: Your Clinical Query Assistant")
st.markdown(
    """
    **DrInfo.ai** is a search engine designed to assist medical professionals and researchers 
    by providing concise, summarized answers to clinical queries. Ask your questions and retrieve 
    answers powered by PubMed and other trusted sources.
    """
)


def submit_query(user_input: str):
    """
    Submits the user's clinical query to the FastAPI backend and updates the session state with the response.

    Args:
        user_input (str): The clinical question entered by the user.

    Raises:
        requests.exceptions.RequestException: If there is a problem connecting to the backend.
    """
    # Clear previous chat history
    st.session_state.messages = [{"role": "user", "content": user_input}]

    # Call FastAPI backend
    with st.spinner("Retrieving information, please wait..."):
        try:
            # Send the query to the backend
            response = requests.post(API_URL, json={"user_message": user_input})
            response.raise_for_status()  # Raise an error for HTTP issues

            # Parse response
            data = response.json()
            st.session_state.messages.append(
                {"role": "bot", "content": data.get("summary", "No summary provided.")}
            )
        except requests.exceptions.RequestException as e:
            st.error(f"Error connecting to the backend: {e}")
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")


def display_chat_history():
    """
    Displays the chat history (user queries and bot responses) on the Streamlit frontend.
    """
    st.markdown("### Search Results")
    for message in st.session_state.messages:
        if message["role"] == "user":
            st.markdown(f"**You:** {message['content']}")
        elif message["role"] == "bot":
            st.markdown(f"**DrInfo.ai:** {message['content']}")


def main():
    """
    Main function for rendering the Streamlit frontend.
    """
    # User Input Section
    with st.container():
        st.markdown("### Ask Your Clinical Question")
        user_input = st.text_input(
            "Enter your question here:", placeholder="e.g., What are the latest treatments for diabetes?"
        )
        submit_button = st.button("Submit")

    # Query Submission Logic
    if submit_button:
        if user_input.strip():
            submit_query(user_input)
        else:
            st.warning("Please enter a valid query!")

    # Display Chat History
    display_chat_history()

    # Footer Section
    st.markdown(
        """
        ---
        **Disclaimer:** This tool provides information from publicly available sources.
        It does not replace professional medical advice. Always consult a healthcare provider for medical concerns.
        """
    )


if __name__ == "__main__":
    main()