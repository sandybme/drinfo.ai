"""
main.py

This module initializes the FastAPI app for the DrInfo.AI project. It includes the setup for serving
static files, routing requests to the PubMed router, and providing utility endpoints such as the favicon.

Modules:
- `pubmed`: Handles PubMed-related queries.
"""

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from app.routers import pubmed
import os

# Initialize FastAPI app
app = FastAPI()

# Include the PubMed router
app.include_router(pubmed.router, prefix="/api", tags=["PubMed"])

# Mount the static directory to serve favicon and other static files
static_dir = os.path.join(os.path.dirname(__file__), "static")
if not os.path.exists(static_dir):
    os.makedirs(static_dir)  # Create the static directory if it doesn't exist
app.mount("/static", StaticFiles(directory=static_dir), name="static")


@app.get("/favicon.ico", include_in_schema=False)
async def favicon():
    """
    Serve the favicon.ico file.

    Returns:
        FileResponse: The favicon.ico file if found.
        dict: A message indicating the favicon is not found.

    Example:
        If the `favicon.ico` file is placed in the `static` directory, 
        it will be served at the `/favicon.ico` endpoint.
    """
    favicon_path = os.path.join(static_dir, "favicon.ico")
    if os.path.exists(favicon_path):
        return FileResponse(favicon_path)
    return {"message": "Favicon not found"}


@app.get("/")
def read_root():
    """
    Root endpoint of the application.

    Returns:
        dict: A welcome message for users accessing the API.

    Example:
        Accessing the root endpoint via `http://localhost:8000/` will return:
        {
            "message": "Welcome to the Doctor's Query API"
        }
    """
    return {"message": "Welcome to the Doctor's Query API"}