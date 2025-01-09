from sqlalchemy import create_engine, Column, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import json

# Create a base class for the ORM models
Base = declarative_base()

# Define the database connection and session
DATABASE_URL = "sqlite:///./test.db"  # This is for local SQLite, change it to your database URL for production
engine = create_engine(DATABASE_URL, connect_args={"check_same_thread": False})
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Define the PubMedArticle model
class PubMedArticle(Base):
    __tablename__ = "pubmed_articles"

    pubmed_id = Column(String, primary_key=True, index=True)
    title = Column(String)
    abstract = Column(String)
    link = Column(String)
    json_data = Column(String)  # Storing the formatted JSON

    def __repr__(self):
        return f"<PubMedArticle(pubmed_id={self.pubmed_id}, title={self.title}, link={self.link})>"

# Create the tables in the database
Base.metadata.create_all(bind=engine)

# Function to insert article into the database
import json
from sqlalchemy.orm import Session
from app.database import PubMedArticle, SessionLocal

def insert_pubmed_article(formatted_json):
    db = SessionLocal()
    try:
        articles = json.loads(formatted_json)  # Parse the formatted JSON string
        
        for article in articles:
            pubmed_id = article.get('pubmed_id')
            title = article.get('title')
            abstract = article.get('abstract')  # This is a list
            link = article.get('link')
            
            # Convert abstract list to string
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
        db.close()