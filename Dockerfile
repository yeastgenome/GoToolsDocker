# =========================================
# GO Tools Docker - FastAPI Python Implementation
# =========================================
# This is a pure Python implementation that replaces the previous
# Perl-based GO-TermFinder with Python using scipy for statistics.
# =========================================

FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
        wget; \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Create app directory structure
RUN set -eux; \
    install -d -m 755 /var/www/app; \
    install -d -m 1777 /var/www/tmp; \
    install -d -m 755 /var/www/data; \
    install -d -m 755 /var/www/cache

WORKDIR /var/www/app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY www/FlaskApp/FlaskApp/*.py ./
COPY www/data/gene_ontology.obo /var/www/data/
COPY www/data/gene_association.sgd /var/www/data/
COPY www/data/slim_*.lst /var/www/data/

# Environment variables
ENV DATA_DIR=/var/www/data/ \
    TMP_DIR=/var/www/tmp/ \
    CACHE_DIR=/var/www/cache/ \
    S3_BUCKET=""

# Health check
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
    CMD curl -fsS http://localhost:8000/ || exit 1

EXPOSE 8000

# Run with uvicorn
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
