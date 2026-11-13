# Multi-stage or clean single-stage lightweight Python image
FROM python:3.11-slim

# Set environment variables for Python in container
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PORT=8080

# Set working directory inside container
WORKDIR /app

# Install compilation tools required by reedsolo C-extension (falls back to pure python if needed)
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy project source files
COPY . .

# Ensure storage directories exist
RUN mkdir -p dashboard_storage/uploads dashboard_storage/exports dashboard_storage/recovered

# Expose standard container port
EXPOSE 8080

# Start HELIX-CORE bi-directional dashboard
CMD ["python", "dashboard.py"]
