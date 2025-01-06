from openai import OpenAI
import os
from dotenv import load_dotenv

load_dotenv()
# Initialize OpenAI client
client = OpenAI()
client.api_key = os.getenv("OPENAI_API_KEY")

def parse_query_with_llm(doctor_query: str) -> str:
    """Convert a natural language query into a structured PubMed query."""
    prompt = f"""
    A doctor is looking for clinical evidence. Convert the following question into a PubMed query:
    Question: "{doctor_query}"
    Provide the PubMed query syntax. Just the syntax alone. Strictly no extra text!
    """
    response = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "You are a helpful assistant for PubMed query formulation."},
            {"role": "user", "content": prompt}
        ]
    )
    return response.choices[0].message.content.strip()

def summarize_with_llm(articles: list) -> str:
    """Summarize findings from PubMed articles."""
    content = "\n".join([
        f"Title: {article['title']}\nAbstract: {article['abstract']}\nLink: {article['link']}"
        for article in articles
    ])
    prompt = f"""
    Summarize the findings from the following PubMed articles. Include references:
    {content}
    """
    response = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "Summarize PubMed findings and include references as links."},
            {"role": "user", "content": prompt}
        ]
    )
    return response.choices[0].message.content.strip()