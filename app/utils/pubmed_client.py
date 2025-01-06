from Bio import Entrez
import os
from dotenv import load_dotenv
load_dotenv()
Entrez.email = os.getenv("EMAIL")

def search_pubmed(query, max_results=10):
    """Search PubMed for articles."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, retmode="xml")
    results = Entrez.read(handle)
    handle.close()
    return results.get('IdList', [])

def fetch_details(pubmed_ids):
    """Fetch detailed metadata for PubMed IDs."""
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

def format_results(articles):
    """Format articles into a structured response."""
    base_url = "https://pubmed.ncbi.nlm.nih.gov/"
    results = []
    for article in articles.get("PubmedArticle", []):
        pubmed_id = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", "No abstract available")
        link = f"{base_url}{pubmed_id}/"
        results.append({"pubmed_id": pubmed_id, "title": title, "abstract": abstract, "link": link})
    return results