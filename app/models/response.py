from pydantic import BaseModel
from typing import List

class Article(BaseModel):
    """
    Represents a PubMed article.

    Attributes:
        title (str): Title of the article.
        abstract (str): Abstract or summary of the article.
        link (str): URL to the article on PubMed.
    """
    title: str
    abstract: str
    link: str

class QueryResponse(BaseModel):
    """
    Represents the response to a clinical query.

    Attributes:
        pubmed_query (str): The structured PubMed query formulated by an LLM.
        summary (str): The summarized findings generated from the query results.
    """
    pubmed_query: str
    summary: str