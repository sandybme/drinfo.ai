from pydantic import BaseModel

class QueryRequest(BaseModel):
    """
    Represents the request structure for querying the API.

    Attributes:
        user_message (str): The clinical question or input message from the user.
    """
    user_message: str  # The clinical query or message input from the user