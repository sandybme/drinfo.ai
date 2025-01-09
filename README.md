# drinfo.ai
DrInfo.ai

DrInfo.ai helps medical professionals and researchers get answers to clinical queries from trusted sources like PubMed. This project provides a web application where users can submit clinical queries and retrieve summarized findings from PubMed articles.

---

## .env.example

Before running the application, you should create a `.env` file based on the `.env.example` provided. This file contains environment variables needed for the application to function correctly. (Rename this file to just .env)

Create a `.env` file in the root directory and update it with your specific details.

```env
# .env.example
API_URL_LOCAL=http://127.0.0.1:8000/api/chat  # Local backend URL
API_URL_PROD=http://your-production-backend-url/api/chat  # Production backend URL
ENV=P  # Set to 'L' for local environment or 'P' for production environment
OPENAI_API_KEY=your_openai_api_key  # Your OpenAI API key
EMAIL=your_email  # Your email for PubMed requests
```

## Running the applicaton

## 1. Using docker in production server

To build and push the Docker image for production, run the following commands:

### Build the Docker Image

```bash
docker buildx build --platform linux/amd64,linux/arm64 -t docker-user-id/drinfo-ai:latest --build-arg ENV=P --push .
```

- `--platform` specifies the platform for which the image is built (amd64, arm64).
- `--build-arg ENV=P` sets the environment to production during the build.
- `--push` pushes the image to Docker Hub after building.

### Pull and Run Docker Image on EC2 Instance or remote platforms

- **Pull the Docker Image** from Docker Hub:
```bash
docker pull docker-user-id/drinfo-ai:latest
```

- **Run the Docker container** in detached mode:
```bash
docker run -d -p 8000:8000 -p 8501:8501 -e ENV=P --name drinfo-ai docker-user-id/drinfo-ai:latest
```
- `-d` runs the container in detached mode.
- `-p` maps the container ports to the host.
- `-e ENV=P` sets the environment to production.



### Stopping the Docker Container
```bash
docker stop drinfo-ai
```

### Removing the Docker Container
```bash
docker rm drinfo-ai
```

## 2. Running the App Locally

### Build the Docker Image

```bash
docker buildx build --platform linux/amd64,linux/arm64 -t docker-user-id/drinfo-ai:latest --build-arg ENV=L --push .
```

- `--platform` specifies the platform for which the image is built (amd64, arm64).
- `--build-arg ENV=L` sets the environment to local during the build.
- `--push` pushes the image to Docker Hub after building.

### Using Docker Container:
- **Run the Docker container** in detached mode:
```bash
docker run -d -p 8000:8000 -p 8501:8501 -e ENV=L --name drinfo-ai docker-user-id/drinfo-ai:latest
```
- `-d` runs the container in detached mode.
- `-p` maps the container ports to the host.
- `-e ENV=L` sets the environment to local.

### If you wish to run the app locally using FastAPI and Streamlit:

- Install the dependencies by running:
```bash
git clone https://github.com/sandybme/drinfo.ai.git
cd drinfo.ai
pip install -r requirements.txt
```

- Start FastAPI:
```bash
uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
```

- Start Streamlit:
```bash
streamlit run app/frontend.py
```

- Open the app in your browser: `http://localhost:8501`



