from fastapi import APIRouter, HTTPException
from app.models.request import QueryRequest
from app.utils.openai_client import parse_query_with_llm

router = APIRouter()

@router.post("/")
def parse_query(request: QueryRequest):
    """Parse a natural language query into a PubMed query."""
    try:
        pubmed_query = parse_query_with_llm(request.user_message)
        return {"pubmed_query": pubmed_query}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))