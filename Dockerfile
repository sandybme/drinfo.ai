# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory inside the container
WORKDIR /app

# Copy only the necessary files to leverage Docker's cache
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code to the container
COPY . .

# Copy the supervisord configuration
COPY supervisord.conf /etc/supervisord.conf

# Expose ports for FastAPI and Streamlit
EXPOSE 8000 8501

# Start both FastAPI and Streamlit using supervisord
CMD ["supervisord", "-c", "/etc/supervisord.conf"]