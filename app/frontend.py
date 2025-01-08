import streamlit as st
import requests
import os
from dotenv import load_dotenv
import re
load_dotenv()

# API URL (Local or Prod)
ENV = os.getenv("ENV")
ENV = "L"
if ENV == "L":
    API_URL = os.getenv("API_URL_LOCAL")
else:
    API_URL = os.getenv("API_URL_PROD")
# print(API_URL)
# Streamlit App Configuration
st.set_page_config(page_title="DrInfo.ai", layout="wide")

# Initialize session state for messages
if "messages" not in st.session_state:
    st.session_state.messages = []

# App Title and Description
st.title("DrInfo.ai: Your Clinical Query Assistant")
st.markdown(
    """
    **DrInfo.ai** helps medical professionals and researchers get answers to clinical queries 
    from trusted sources like PubMed. Ask a question, and we will summarize the findings for you.
    """
)

# User Input Section
with st.container():
    st.markdown("### Ask Your Clinical Question")
    user_input = st.text_input("Enter your question here:", placeholder="e.g., What are the latest treatments for diabetes?")
    submit_button = st.button("Search")

# Query Submission Logic
if submit_button:
    if user_input.strip():
        # Step 1: Display User Query in the Chat
        st.session_state.messages = [{"role": "user", "content": user_input}]
        # if "full_summary" not in st.session_state:
        st.session_state.full_summary = ""
        
        # Step 2: Set up placeholders for each step
        pubmed_query_placeholder = st.empty()
        retrieving_placeholder = st.empty()
        summary_placeholder = st.empty()

        # Step 3: Request backend to process the query
        with st.spinner("Processing..."):
            try:
                response = requests.post(f"{API_URL}", json={"user_message": user_input}, stream=True)
                
                for line in response.iter_lines():
                    if line:
                        message = line.decode("utf-8").strip()
                        message = re.sub(r'^data[:\s]*', '', message).strip()
                        if message.startswith("Formulating PubMed query"):
                            message = message.replace('Formulating PubMed query:', '').strip() 
                            pubmed_query_placeholder.markdown(f"**Formulated PubMed query:** {message}")
                        elif message.startswith("Searching PubMed"):
                            retrieving_placeholder.markdown(f"**Retrieving articles:** {message}")
                        elif message:  # This checks if the message is not empty
                        # Append the current message to the full summary in session state
                            st.session_state.full_summary += message + "\n"
                            # Update the summary placeholder with the full summary
                            summary_placeholder.markdown(f"**Summary:**\n{st.session_state.full_summary}")
                            print(message)
            except requests.exceptions.RequestException as e:
                st.error(f"Error connecting to the backend: {e}")
            except Exception as e:
                st.error(f"An unexpected error occurred: {e}")
    else:
        st.warning("Please enter a valid query!")

# Footer Section
st.markdown(
    """
    ---
    **Disclaimer:** This tool provides information from publicly available sources. It does not replace professional medical advice. Always consult a healthcare provider for medical concerns.
    """
)