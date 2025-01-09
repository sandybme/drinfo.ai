from fastapi import APIRouter, HTTPException
from sse_starlette.sse import EventSourceResponse
from app.utils.pubmed_client import search_pubmed, fetch_details, format_results
from app.utils.openai_client import parse_query_with_llm, summarize_with_llm
from app.models.request import QueryRequest
from app.models.response import QueryResponse
from app.database import insert_pubmed_article  # Importing insert function from database.py
import time
import logging
import json

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

router = APIRouter()

async def generate_stream(request: QueryRequest):
    """
    Asynchronous generator to simulate data streaming for PubMed query formulation and summarization,
    while logging the time taken for each step.
    """
    try:
        # Step 1: Formulate PubMed query
        start_time = time.time()  # Start time for Step 1
        pubmed_query = parse_query_with_llm(request.user_message)
        pubmed_query_time = time.time() - start_time  # Time taken for Step 1
        logger.info(f"Step 1 completed in {pubmed_query_time:.2f} seconds.")
        yield f"Formulating PubMed query: {pubmed_query}\n\n"

        # Step 2: Search PubMed
        start_time = time.time()  # Start time for Step 2
        pubmed_ids = search_pubmed(pubmed_query)
        search_time = time.time() - start_time  # Time taken for Step 2

        if not pubmed_ids:
            yield "No relevant articles found.\n\n"
            logger.info(f"Step 2 completed in {search_time:.2f} seconds. No articles found.")
            return

        logger.info(f"Step 2 completed in {search_time:.2f} seconds. Found {len(pubmed_ids)} articles.")
        yield "Searching PubMed database...\n\n"

        # Step 3: Fetch PubMed details
        start_time = time.time()  # Start time for Step 3
        articles = fetch_details(pubmed_ids)
        formatted_articles = format_results(articles)

        # Send articles to the database (Database handling happens in database.py)
        formatted_json = json.dumps(formatted_articles, indent=2)
        insert_pubmed_article(formatted_json)  # Inserting articles into the database

        fetch_time = time.time() - start_time  # Time taken for Step 3
        logger.info(f"Step 3 completed in {fetch_time:.2f} seconds. Fetched {len(articles)} articles.")
        yield "Fetching PubMed details...\n\n"

        # Step 4: Summarize results
        start_time = time.time()  # Start time for Step 4
        summary = summarize_with_llm(formatted_articles, request.user_message)
        summary_time = time.time() - start_time  # Time taken for Step 4

        logger.info(f"Step 4 completed in {summary_time:.2f} seconds. Summary generated.")
        yield f"{summary}\n\n"

    except Exception as e:
        logger.error(f"Error during stream processing: {e}")
        yield f"Error occurred: {str(e)}\n\n"

@router.post("/chat", response_model=QueryResponse)
async def chat_with_doctor(request: QueryRequest):
    """
    Endpoint to process clinical queries, stream the process steps, and return summarized results.
    """
    try:
        return EventSourceResponse(generate_stream(request))  # No need to pass db here
    except Exception as e:
        logger.error(f"Error in processing request: {e}")
        raise HTTPException(status_code=500, detail="Error occurred while processing the query.")