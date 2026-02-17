# GO Tools API

A FastAPI-based REST API for Gene Ontology (GO) analysis tools:
- **GO Term Finder**: Find enriched GO terms in a gene list using hypergeometric test
- **GO Slim Mapper**: Map genes to broader GO Slim categories

## Features

- Pure Python implementation (no Perl dependencies)
- Hypergeometric statistical test for enrichment analysis
- Multiple testing correction (Bonferroni, Benjamini-Hochberg FDR)
- FastAPI with automatic OpenAPI/Swagger documentation
- Compatible with SGD (Saccharomyces Genome Database) gene annotations

---

## Option 1: Running with Docker (Recommended)

### Prerequisites
- Docker and Docker Compose installed

### Quick Start

```bash
# Clone the repository
git clone https://github.com/yeastgenome/GoToolsDocker.git
cd GoToolsDocker

# Build and run
docker-compose up --build

# Or build and run manually
docker build -t go-tools .
docker run -p 8000:8000 go-tools
```

### With S3 Upload Support

```bash
docker run -p 8000:8000 \
  -e S3_BUCKET=your-bucket-name \
  -e AWS_ACCESS_KEY_ID=your-key \
  -e AWS_SECRET_ACCESS_KEY=your-secret \
  go-tools
```

### Access the API

- API: http://localhost:8000
- Swagger UI: http://localhost:8000/docs
- OpenAPI JSON: http://localhost:8000/openapi.json

---

## Option 2: Running on EC2 (Without Docker)

### Prerequisites

- Python 3.11+
- pip

### Step 1: Install System Dependencies

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install -y python3 python3-pip python3-venv curl wget

# Amazon Linux 2023
sudo dnf install -y python3.11 python3.11-pip curl wget
```

### Step 2: Clone and Setup

```bash
# Clone repository
git clone https://github.com/yeastgenome/GoToolsDocker.git
cd GoToolsDocker

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Step 3: Setup Data Directory

```bash
# Create directories
sudo mkdir -p /var/www/data /var/www/tmp /var/www/cache
sudo chown -R $USER:$USER /var/www

# Copy data files
cp www/data/gene_ontology.obo /var/www/data/
cp www/data/gene_association.sgd /var/www/data/
cp www/data/slim_*.lst /var/www/data/
```

### Step 4: Run the Application

```bash
# Set environment variables
export DATA_DIR=/var/www/data/
export TMP_DIR=/var/www/tmp/
export CACHE_DIR=/var/www/cache/
export S3_BUCKET=""  # Optional: set your S3 bucket name

# Run with uvicorn
cd www/FlaskApp/FlaskApp
uvicorn main:app --host 0.0.0.0 --port 8000
```

### Step 5: Run as a Service (Optional)

Create a systemd service file:

```bash
sudo tee /etc/systemd/system/go-tools.service > /dev/null <<EOF
[Unit]
Description=GO Tools API
After=network.target

[Service]
User=$USER
WorkingDirectory=/home/$USER/GoToolsDocker/www/FlaskApp/FlaskApp
Environment="DATA_DIR=/var/www/data/"
Environment="TMP_DIR=/var/www/tmp/"
Environment="CACHE_DIR=/var/www/cache/"
Environment="S3_BUCKET="
ExecStart=/home/$USER/GoToolsDocker/venv/bin/uvicorn main:app --host 0.0.0.0 --port 8000
Restart=always

[Install]
WantedBy=multi-user.target
EOF

# Enable and start the service
sudo systemctl daemon-reload
sudo systemctl enable go-tools
sudo systemctl start go-tools

# Check status
sudo systemctl status go-tools
```

### Step 6: Setup Nginx Reverse Proxy (Optional)

```bash
sudo apt install -y nginx

sudo tee /etc/nginx/sites-available/go-tools > /dev/null <<EOF
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
    }
}
EOF

sudo ln -s /etc/nginx/sites-available/go-tools /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

---

## API Endpoints

### Health Check
```
GET /
```

### GO Term Finder (Full)
```
GET/POST /gotermfinder
```
Parameters:
- `genes` (required): Pipe-separated gene names (e.g., `CYT1|QCR2|ACT1`)
- `aspect`: GO aspect - `P` (Process), `F` (Function), `C` (Component). Default: `F`
- `pvalue`: P-value cutoff. Default: `0.01`
- `genes4bg`: Pipe-separated background genes (optional)
- `evidence`: Pipe-separated evidence codes (optional)
- `FDR`: Enable FDR correction (optional)

Example:
```bash
curl "http://localhost:8000/gotermfinder?genes=CYT1|QCR2|ACT1&aspect=F&pvalue=0.01"
```

### GO Term Finder (Simple JSON)
```
GET/POST /termfinder
```
Returns a simple JSON array of enriched terms.

Example:
```bash
curl "http://localhost:8000/termfinder?genes=CYT1|QCR2|ACT1&aspect=F"
```

### GO Slim Mapper
```
GET/POST /goslimmapper
```
Parameters:
- `genes` (required): Pipe-separated gene names
- `aspect` (required): GO aspect - `P`, `F`, or `C`
- `terms` (required): Pipe-separated slim terms (format: `term name ; GO:XXXXXXX`)

Example:
```bash
curl "http://localhost:8000/goslimmapper?genes=CYT1|QCR2|ACT1&aspect=F&terms=molecular%20function;GO:0003674|catalytic%20activity;GO:0003824"
```

---

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `DATA_DIR` | Directory containing OBO and GAF files | `/var/www/data/` |
| `TMP_DIR` | Directory for temporary output files | `/var/www/tmp/` |
| `CACHE_DIR` | Directory for ontology cache | `/var/www/cache/` |
| `S3_BUCKET` | S3 bucket for file uploads (optional) | `` |

---

## Data Files

The following data files are required in `DATA_DIR`:

- `gene_ontology.obo` - Gene Ontology in OBO format
- `gene_association.sgd` - Gene associations in GAF format
- `slim_process.lst` - GO Slim terms for Biological Process
- `slim_function.lst` - GO Slim terms for Molecular Function
- `slim_component.lst` - GO Slim terms for Cellular Component

Download latest files from:
- GO Ontology: http://purl.obolibrary.org/obo/go.obo
- SGD Annotations: https://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz

---

## License

MIT License
