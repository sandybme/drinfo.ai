from fastapi import APIRouter, HTTPException
from app.models.request import QueryRequest
from app.utils.openai_client import parse_query_with_llm

# Initialize the router
router = APIRouter()

@router.post("/")
def parse_query(request: QueryRequest):
    """
    Endpoint to convert a natural language query into a structured PubMed query.

    Args:
        request (QueryRequest): The user-provided natural language query.

    Returns:
        dict: A dictionary containing the parsed PubMed query.

    Raises:
        HTTPException: If any error occurs during query parsing.

    Workflow:
        1. Receive a natural language query from the user.
        2. Use an LLM to parse the query into PubMed-compatible syntax.
        3. Return the structured PubMed query.
    """
    try:
        pubmed_query = parse_query_with_llm(request.user_message)
        return {"pubmed_query": pubmed_query}
    except Exception as e:
        # Raise an HTTP 500 error with the exception details
        raise HTTPException(status_code=500, detail=f"Error parsing query: {str(e)}")