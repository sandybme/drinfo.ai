from Bio import Entrez
import os
from dotenv import load_dotenv
from langchain_core.vectorstores import InMemoryVectorStore
from langchain_openai import OpenAIEmbeddings
import json

# Load environment variables
load_dotenv()
Entrez.email = os.getenv("EMAIL")

class Document:
    """Represents a document with an abstract, metadata, and a unique identifier."""
    
    def __init__(self, page_content: str, metadata: dict, doc_id: str):
        self.page_content = page_content  # The abstract text
        self.metadata = metadata  # A dictionary with link, pubmed_id, and title
        self.id = doc_id  # Unique document identifier
    
    def __repr__(self):
        """String representation of the Document object."""
        return f"Document(id={self.id}, page_content={self.page_content[:100]}..., metadata={self.metadata})"

def search_pubmed(query: str, max_results: int = 10) -> list:
    """
    Search PubMed for article IDs based on a query.
    
    Args:
        query (str): Search term or query.
        max_results (int): Maximum number of results to retrieve.

    Returns:
        list: List of PubMed IDs.
    """
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, retmode="xml")
    results = Entrez.read(handle)
    handle.close()
    return results.get('IdList', [])

def fetch_details(pubmed_ids: list) -> dict:
    """
    Fetch metadata for PubMed articles by IDs.
    
    Args:
        pubmed_ids (list): List of PubMed IDs.

    Returns:
        dict: Detailed metadata for the articles.
    """
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

def format_results(articles: dict) -> list:
    """
    Format PubMed article metadata into structured responses.
    
    Args:
        articles (dict): PubMed article metadata.

    Returns:
        list: List of formatted article details with title, abstract, and link.
    """
    base_url = "https://pubmed.ncbi.nlm.nih.gov/"
    results = []
    for article in articles.get("PubmedArticle", []):
        pubmed_id = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", "No abstract available")
        link = f"{base_url}{pubmed_id}/"
        results.append({
            "pubmed_id": pubmed_id,
            "title": title,
            "abstract": abstract,
            "link": link
        })
    return results

def rerank_results(articles: list, user_query: str) -> list:
    """
    Rerank PubMed articles based on the user query using embeddings and cosine similarity.
    
    Args:
        articles (list): List of formatted article metadata.
        user_query (str): The query provided by the user.

    Returns:
        list: Sorted PubMed article IDs based on relevance to the query.
    """
    embeddings = OpenAIEmbeddings(model="text-embedding-3-large")
    
    documents = [
        Document(
            page_content=" ".join(article["abstract"]),  # Abstract text as page content
            metadata={
                'source': article['link'],  # Source link
                'pubmed_id': article['pubmed_id'],  # PubMed ID
                'title': article['title']  # Article Title
            },
            doc_id=article['pubmed_id']  # PubMed ID as document ID
        )
        for article in articles
    ]

    # Create vector store from documents and query embeddings
    vectorstore = InMemoryVectorStore.from_documents(documents=documents, embedding=embeddings)

    # Use the vector store as a retriever and retrieve top-k documents
    retriever = vectorstore.as_retriever(search_kwargs={'k': 20})
    retrieved_documents = retriever.invoke(user_query)

    # Extract sorted PubMed IDs from the retrieved documents
    reranked_ids = [doc.id for doc in retrieved_documents]
    
    return reranked_ids