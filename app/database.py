from sqlalchemy import create_engine, Column, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import json
from sqlalchemy.orm import Session

# Create a base class for the ORM models
Base = declarative_base()

# Define the database connection and session
DATABASE_URL = "sqlite:///./test.db"  # Local SQLite for testing; change it to your database URL for production
engine = create_engine(DATABASE_URL, connect_args={"check_same_thread": False})
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Define the PubMedArticle model
class PubMedArticle(Base):
    """Database model representing a PubMed article."""
    
    __tablename__ = "pubmed_articles"

    pubmed_id = Column(String, primary_key=True, index=True)
    title = Column(String)
    abstract = Column(String)
    link = Column(String)
    json_data = Column(String)  # Storing the formatted JSON

    def __repr__(self):
        """String representation of the PubMedArticle."""
        return f"<PubMedArticle(pubmed_id={self.pubmed_id}, title={self.title}, link={self.link})>"

# Create the tables in the database
Base.metadata.create_all(bind=engine)

def insert_pubmed_article(formatted_json: str):
    """
    Insert articles into the database from the provided JSON formatted string.

    Args:
        formatted_json (str): A JSON string containing articles data with fields
                              like pubmed_id, title, abstract, and link.

    Returns:
        None
    """
    db = SessionLocal()
    try:
        # Parse the formatted JSON string
        articles = json.loads(formatted_json)

        for article in articles:
            pubmed_id = article.get('pubmed_id')
            title = article.get('title')
            abstract = article.get('abstract')  # Abstract is expected to be a list
            link = article.get('link')

            # Convert abstract list to string if it's a list
            abstract_str = " ".join(abstract) if isinstance(abstract, list) else abstract

            # Check if article already exists in the database
            existing_article = db.query(PubMedArticle).filter(PubMedArticle.pubmed_id == pubmed_id).first()

            if not existing_article:  # Insert only if the pubmed_id is not already present
                db_article = PubMedArticle(
                    pubmed_id=pubmed_id,
                    title=title,
                    abstract=abstract_str,  # Save abstract as a string
                    link=link,
                    json_data=formatted_json,
                )
                db.add(db_article)
                db.commit()
                db.refresh(db_article)  # Refresh to get the updated values
    finally:
        db.close()  # Ensure the session is closed