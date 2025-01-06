import streamlit as st
import requests

# Backend API URL
API_URL = "http://127.0.0.1:8000/api/chat"

# Streamlit App Configuration
st.set_page_config(page_title="DrInfo.ai", layout="wide")

# Initialize session state for messages
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

# User Input Section
with st.container():
    st.markdown("### Ask Your Clinical Question")
    user_input = st.text_input("Enter your question here:", placeholder="e.g., What are the latest treatments for diabetes?")
    submit_button = st.button("Submit")

# Query Submission Logic
if submit_button:
    if user_input.strip():
        # Clear previous chat history
        st.session_state.messages = [{"role": "user", "content": user_input}]

        # Call FastAPI backend
        with st.spinner("Retrieving information, please wait..."):
            try:
                response = requests.post(API_URL, json={"user_message": user_input})
                response.raise_for_status()

                data = response.json()
                st.session_state.messages.append({"role": "bot", "content": data.get("summary", "No summary provided.")})
            except requests.exceptions.RequestException as e:
                st.error(f"Error connecting to the backend: {e}")
            except Exception as e:
                st.error(f"An unexpected error occurred: {e}")
    else:
        st.warning("Please enter a valid query!")

# Display Chat History
st.markdown("### Search Results")
for message in st.session_state.messages:
    if message["role"] == "user":
        st.markdown(f"**You:** {message['content']}")
    elif message["role"] == "bot":
        st.markdown(f"**DrInfo.ai:** {message['content']}")

# Footer Section
st.markdown(
    """
    ---
    **Disclaimer:** This tool provides information from publicly available sources. It does not replace professional medical advice. Always consult a healthcare provider for medical concerns.
    """
)   