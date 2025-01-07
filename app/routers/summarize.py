from fastapi import APIRouter, HTTPException
from app.utils.openai_client import summarize_with_llm

# Initialize the router
router = APIRouter()

@router.post("/")
def summarize_articles(articles: list):
    """
    Endpoint to generate a summary from PubMed articles.

    Args:
        articles (list): A list of dictionaries containing details of PubMed articles.
                         Each dictionary should include keys like 'title', 'abstract', and 'link'.

    Returns:
        dict: A dictionary containing the summarized findings from the provided articles.

    Raises:
        HTTPException: If any error occurs during summarization.

    Workflow:
        1. Receive a list of articles as input.
        2. Use an LLM to generate a concise summary of the findings.
        3. Return the generated summary as a JSON response.
    """
    try:
        # Generate summary using LLM
        summary = summarize_with_llm(articles)
        return {"summary": summary}
    except Exception as e:
        # Log and raise HTTP 500 error with exception details
        raise HTTPException(status_code=500, detail=f"Error summarizing articles: {str(e)}")