import logging
from fastapi import APIRouter, HTTPException
from app.utils.pubmed_client import search_pubmed, fetch_details, format_results
from app.utils.openai_client import parse_query_with_llm, summarize_with_llm
from app.models.request import QueryRequest
from app.models.response import QueryResponse

# Create a logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

router = APIRouter()

@router.post("/chat", response_model=QueryResponse)
def chat_with_doctor(request: QueryRequest):
    """
    Chat endpoint for doctors to ask questions and retrieve PubMed results.
    """
    try:
        logger.info("Received request: %s", request.user_message)   

        # Step 1: Formulate PubMed query
        pubmed_query = parse_query_with_llm(request.user_message)
        # logger.info("Generated PubMed query: %s", pubmed_query)

        # Step 2: Search PubMed
        pubmed_ids = search_pubmed(pubmed_query)
        # print(pubmed_ids)
        if not pubmed_ids:
            # logger.warning("No PubMed articles found for query: %s", pubmed_query)
            return QueryResponse(bot_response="No relevant articles found.", results=[])
        # logger.info("Retrieved PubMed IDs: %s", pub   med_ids)

        # Step 3: Fetch PubMed details
        articles = fetch_details(pubmed_ids)
        formatted_articles = format_results(articles)
        # print(articles)
        # logger.info("Fetched PubMed articles: %s", articles)

        # Step 4: Summarize results
        summary = summarize_with_llm(formatted_articles)
        logger.info("Generated summary: %s", summary)

        # Return response
        return QueryResponse(pubmed_query=pubmed_query, summary=summary)

    except Exception as e:
        logger.error("An error occurred: %s", str(e), exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))