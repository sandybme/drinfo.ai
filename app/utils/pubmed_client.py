from Bio import Entrez
import os
from dotenv import load_dotenv
from langchain_core.vectorstores import InMemoryVectorStore
from langchain_openai import OpenAIEmbeddings


# Load environment variables
load_dotenv()
Entrez.email = os.getenv("EMAIL")

class Document:
    def __init__(self, page_content: str, metadata: dict, doc_id: str):
        self.page_content = page_content  # The abstract text
        self.metadata = metadata  # A dictionary with link, pubmed_id, and title
        self.id = doc_id  # Adding a unique ID field
    
    def __repr__(self):
        return f"Document(id={self.id}, page_content={self.page_content[:100]}..., metadata={self.metadata})"  # Show first 100 chars of page_content


def search_pubmed(query, max_results=10):
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

def fetch_details(pubmed_ids):
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

def format_results(articles):
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

def rerank_results(articles,user_query):

    embeddings = OpenAIEmbeddings(model="text-embedding-3-large")
    
    documents = [
    Document(
        page_content=" ".join(article["abstract"]),  # The abstract as page content
        metadata={
            'source': article['link'],  # The link as source
            'pubmed_id': article['pubmed_id'],  # PubMed ID
            'title': article['title']  # Article Title
        },
        doc_id= article['pubmed_id']  # Using PubMed ID as the unique ID for the document
    )

    for article in articles ]
    # print(len(documents))

    vectorstore = InMemoryVectorStore.from_documents(
    documents=documents,
    embedding=embeddings,
)
# # Use the vectorstore as a retriever
    retriever = vectorstore.as_retriever(search_kwargs={'k': 20})
    retrieved_documents = retriever.invoke(user_query)
    # print(len(retrieved_documents))
    reranked_ids = [article.id for article in retrieved_documents]
    return reranked_ids
