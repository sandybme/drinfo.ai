from openai import OpenAI
from Bio import Entrez
from dotenv import load_dotenv
import os
import logging

# --- Setup Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# --- Load Environment Variables ---
load_dotenv()

# Initialize OpenAI client
client = OpenAI()
openai_api_key = os.getenv("OPENAI_API_KEY")
email = os.getenv("EMAIL")

# Configure OpenAI and Entrez
client.api_key = openai_api_key
Entrez.email = email

PUBMED_BASE_URL = "https://pubmed.ncbi.nlm.nih.gov/"


# --- Utility Functions ---
def handle_openai_response(response):
    """Safely extract content from OpenAI's response."""
    try:
        return response.choices[0].message.content.strip()
    except (KeyError, IndexError) as e:
        logging.error(f"Failed to parse OpenAI response: {e}")
        raise RuntimeError("Invalid response from OpenAI API.") from e


# --- Step 1: Parse Query Using OpenAI ---
def parse_query_with_llm(doctor_query):
    """Convert a natural language query into a structured PubMed query."""
    logging.info("Parsing query with OpenAI LLM...")
    prompt = f"""
    A doctor is looking for clinical evidence. Convert the following question into a PubMed query:
    Question: "{doctor_query}"
    Provide the PubMed query syntax.
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "You are a helpful assistant for PubMed query formulation."},
                {"role": "user", "content": prompt}
            ]
        )
        return handle_openai_response(response)
    except Exception as e:
        logging.error(f"Error while querying OpenAI API: {e}")
        raise


# --- Step 2: Fetch Data from PubMed ---
def search_pubmed(query, max_results=10):
    """Search PubMed for articles matching the query."""
    logging.info("Searching PubMed...")
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, retmode="xml")
        results = Entrez.read(handle)
        handle.close()
        return results.get('IdList', [])
    except Exception as e:
        logging.error(f"Error while querying PubMed: {e}")
        raise


def fetch_details(pubmed_ids):
    """Fetch detailed metadata for a list of PubMed IDs."""
    logging.info("Fetching details for PubMed IDs...")
    try:
        handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        return records
    except Exception as e:
        logging.error(f"Error while fetching PubMed details: {e}")
        raise


# --- Step 3: Format Results ---
def format_results(articles):
    """Extract titles, abstracts, and generate links for each article."""
    logging.info("Formatting PubMed results...")
    results = []
    for article in articles.get("PubmedArticle", []):
        try:
            pubmed_id = article["MedlineCitation"]["PMID"]
            title = article["MedlineCitation"]["Article"]["ArticleTitle"]
            abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", "No abstract available")
            link = f"{PUBMED_BASE_URL}{pubmed_id}/"
            results.append({"pubmed_id": pubmed_id, "title": title, "abstract": abstract, "link": link})
        except KeyError as e:
            logging.warning(f"Skipping an article due to missing data: {e}")
            continue
    return results


# --- Step 4: Summarize Results with OpenAI ---
def summarize_with_llm(articles):
    """Summarize findings from PubMed articles using OpenAI API."""
    logging.info("Summarizing articles with OpenAI LLM...")
    content = "\n".join([
        f"Title: {article['title']}\nAbstract: {article['abstract']}\nLink: {article['link']}"
        for article in articles
    ])
    prompt = f"""
    Summarize the findings from the following PubMed articles. Include references:
    {content}
    """
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "Summarize PubMed findings and include references as links."},
                {"role": "user", "content": prompt}
            ]
        )
        return handle_openai_response(response)
    except Exception as e:
        logging.error(f"Error while summarizing with OpenAI API: {e}")
        raise


# --- Main Program ---
def main():
    # Step 1: Doctor's Query
    doctor_query = "Which medication is better for diabetes: Metformin or Insulin?"
    logging.info(f"Doctor's Query: {doctor_query}")

    try:
        # Step 2: Parse Query with LLM
        pubmed_query = parse_query_with_llm(doctor_query)
        logging.info(f"Generated PubMed Query: {pubmed_query}")

        # Step 3: Search PubMed
        pubmed_ids = search_pubmed(pubmed_query, max_results=5)
        if not pubmed_ids:
            logging.warning("No articles found.")
            return
        logging.info(f"Found {len(pubmed_ids)} articles.")

        # Step 4: Fetch and Format Results
        articles = fetch_details(pubmed_ids)
        formatted_articles = format_results(articles)

        # Step 5: Summarize Results
        summary = summarize_with_llm(formatted_articles)
        logging.info("Summary of Findings:")
        print(summary)
    except Exception as e:
        logging.error(f"An error occurred: {e}")


# Run the Program
if __name__ == "__main__":
    main()