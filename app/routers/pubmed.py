import logging
from fastapi import APIRouter, HTTPException
from app.utils.pubmed_client import search_pubmed, fetch_details, format_results
from app.utils.openai_client import parse_query_with_llm, summarize_with_llm
from app.models.request import QueryRequest
from app.models.response import QueryResponse

# Configure the logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Initialize the router
router = APIRouter()

@router.post("/chat", response_model=QueryResponse)
def chat_with_doctor(request: QueryRequest):
    """
    Endpoint for processing clinical queries and returning summarized PubMed results.

    Args:
        request (QueryRequest): The user-provided clinical query.

    Returns:
        QueryResponse: A response containing the PubMed query and summarized findings.

    Workflow:
        1. Parse the user query to formulate a PubMed query using an LLM.
        2. Search PubMed for relevant articles.
        3. Fetch and format metadata of retrieved articles.
        4. Generate a summary using an LLM.
    """
    try:
        logger.info("Received request: %s", request.user_message)

        # Formulate PubMed query
        pubmed_query = parse_query_with_llm(request.user_message)

        # Search PubMed
        pubmed_ids = search_pubmed(pubmed_query)
        if not pubmed_ids:
            logger.warning("No PubMed articles found for query: %s", pubmed_query)
            return QueryResponse(pubmed_query=pubmed_query, summary="No relevant articles found.")

        # Fetch PubMed article details
        articles = fetch_details(pubmed_ids)
        formatted_articles = format_results(articles)

        # Summarize results
        summary = summarize_with_llm(formatted_articles)
        logger.info("Generated summary: %s", summary)

        return QueryResponse(pubmed_query=pubmed_query, summary=summary)

    except Exception as e:
        logger.error("An error occurred: %s", str(e), exc_info=True)
        raise HTTPException(status_code=500, detail="An error occurred while processing the query.")