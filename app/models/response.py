from pydantic import BaseModel
from typing import List

class Article(BaseModel):
    title: str
    abstract: str
    link: str

class QueryResponse(BaseModel):
    pubmed_query: str  # The formulated PubMed query
    summary: str  # A summarized response from OpenAI or similar LLM