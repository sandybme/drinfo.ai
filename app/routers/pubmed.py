from fastapi import APIRouter, HTTPException
from sse_starlette.sse import EventSourceResponse
from app.utils.pubmed_client import search_pubmed, fetch_details, format_results
from app.utils.openai_client import parse_query_with_llm, summarize_with_llm
from app.models.request import QueryRequest
from app.models.response import QueryResponse
import time
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

router = APIRouter()

async def generate_stream(request: QueryRequest):
    """
    Asynchronous generator to simulate data streaming for PubMed query formulation and summarization.
    """
    try:
        # Step 1: Formulate PubMed query
        pubmed_query = parse_query_with_llm(request.user_message)
        print(pubmed_query)
        yield f"Formulating PubMed query: {pubmed_query}\n\n"
        
        # Step 2: Search PubMed
        pubmed_ids = search_pubmed(pubmed_query)
        if not pubmed_ids:
            yield "No relevant articles found.\n\n"
            return

        # Step 3: Fetch PubMed details
        articles = fetch_details(pubmed_ids)
        formatted_articles = format_results(articles)
        yield f"Searching PubMed database...\n\n"

        # Step 4: Summarize results
        summary = summarize_with_llm(formatted_articles,request.user_message)
        # print(summary)
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
        return EventSourceResponse(generate_stream(request))
    except Exception as e:
        logger.error(f"Error in processing request: {e}")
        raise HTTPException(status_code=500, detail="Error occurred while processing the query.")