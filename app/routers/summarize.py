from fastapi import APIRouter, HTTPException
from app.utils.openai_client import summarize_with_llm

router = APIRouter()

@router.post("/")
def summarize_articles(articles: list):
    """Summarize findings from PubMed articles."""
    try:
        summary = summarize_with_llm(articles)
        return {"summary": summary}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))