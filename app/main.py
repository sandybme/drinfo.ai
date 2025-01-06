from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from app.routers import pubmed
import os

app = FastAPI()

# Include the PubMed router
app.include_router(pubmed.router, prefix="/api", tags=["PubMed"])

# Mount the static directory to serve favicon and other static files
static_dir = os.path.join(os.path.dirname(__file__), "static")
if not os.path.exists(static_dir):
    os.makedirs(static_dir)  # Create the static directory if it doesn't exist
app.mount("/static", StaticFiles(directory=static_dir), name="static")

# Serve favicon
@app.get("/favicon.ico", include_in_schema=False)
async def favicon():
    """
    Serve the favicon.ico file.
    """
    favicon_path = os.path.join(static_dir, "favicon.ico")
    if os.path.exists(favicon_path):
        return FileResponse(favicon_path)
    return {"message": "Favicon not found"}

@app.get("/")
def read_root():
    """
    Root endpoint.c
    """
    return {"message": "Welcome to the Doctor's Query API"}