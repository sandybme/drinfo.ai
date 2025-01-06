from pydantic import BaseModel

class QueryRequest(BaseModel):
    user_message: str  # The expected field name