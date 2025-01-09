from openai import OpenAI
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Initialize OpenAI client
client = OpenAI()
client.api_key = os.getenv("OPENAI_API_KEY")

def parse_query_with_llm(doctor_query: str) -> str:
    """
    Convert a natural language query into a structured PubMed query.

    Args:
        doctor_query (str): Doctor's clinical question in natural language.

    Returns:
        str: Structured PubMed query syntax or an "" if not a clinical question.
    """
    prompt = f"""
    A doctor is looking for clinical evidence. Convert the following question into a PubMed query:
    Question: "{doctor_query}"
    Provide the PubMed query syntax. Just the syntax alone. Strictly no extra text! 
    If the question is unrelated to clinical queries, do not provide anything and strictly return an "".
    If the query does not relate to treatments, medical conditions, or clinical evidence, leave it empty.

    Do not return the question verbatim, return only the structured PubMed query if it's clinical.
    """
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        temperature= 0.2,
        messages=[
            {"role": "system", "content": "You are a helpful assistant for PubMed query formulation."},
            {"role": "user", "content": prompt}
        ]
    )

    # Retrieve the response and clean it
    pubmed_query = response.choices[0].message.content.strip()

    # If it's empty or a non-clinical question, return an empty string
    if not pubmed_query or " " not in pubmed_query:
        return ""  # Return empty string if invalid or non-clinical query
    print(pubmed_query)
    # If a valid query is formed, return it
    return pubmed_query

def summarize_with_llm(articles: list, user_query: str) -> str:
    """
    Summarize findings from a list of PubMed articles based on a user's query.

    Args:
        articles (list): List of dictionaries containing article metadata.
        user_query (str): The user's clinical query to provide context for the summary.

    Returns:
        str: Summary of findings, including relevant references.
    """
    # Prepare article information and keep the links for citation
    article_links = []
    content = ""
    for i, article in enumerate(articles):
        article_links.append(f"[PubMed Link {i+1}]({article['link']})")  # Create clickable reference
        content += f"Title: {article['title']}\nAbstract: {article['abstract']}\nLink: {article['link']}\n\n"
    
    # Construct the prompt with the user's query for better context
    prompt = f"""
    Summarize the findings from the PubMed articles related to this {user_query}. Exclude articles irrelavant to the query. Include references wherever relevant:
    {content}
    """
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        temperature= 0.2,
        messages=[
            {"role": "system", "content": "Summarize PubMed findings in an elaborate systematic way. Always include references as links."},
            {"role": "user", "content": prompt}
        ]
    )

    # Get the generated summary
    summary = response.choices[0].message.content.strip()

    # # Replace placeholder citations with actual clickable links
    # for i, link in enumerate(article_links):
    #     summary = summary.replace(f"[PubMed Link {i+1}]", link)  # Replace placeholder with actual link

    return summary