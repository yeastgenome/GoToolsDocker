# =========================================
# Stage 1: build GO-TermFinder (builder)
# =========================================
FROM ubuntu:24.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive

# Make apt more resilient + smaller indexes
RUN set -eux; \
  echo 'Acquire::Retries "5";' > /etc/apt/apt.conf.d/80-retries; \
  echo 'Acquire::http::No-Cache "true";' > /etc/apt/apt.conf.d/99no-cache; \
  echo 'Acquire::Languages "none";' > /etc/apt/apt.conf.d/99no-languages

# Build deps (kept only in builder) — split into small installs
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends build-essential gcc g++ make; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends git wget ca-certificates; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends perl libcgi-pm-perl libgd-perl libgraphviz-perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

WORKDIR /tmp/build

# Fetch repo + GO-TermFinder tarball
RUN set -eux; \
  git clone --depth=1 https://github.com/yeastgenome/GoToolsDocker.git /tmp/GoToolsDocker; \
  wget -O GO-TermFinder-0.86.tar.gz https://cpan.metacpan.org/authors/id/S/SH/SHERLOCK/GO-TermFinder-0.86.tar.gz; \
  tar xzf GO-TermFinder-0.86.tar.gz

# Build/install GO-TermFinder into /opt/gtf (only this gets copied to runtime)
WORKDIR /tmp/build/GO-TermFinder-0.86
RUN set -eux; \
  perl Makefile.PL INSTALL_BASE=/opt/gtf; \
  make; \
  make test || true; \
  make install; \
  find /opt/gtf -type f -name '*.pod' -delete || true

# =========================================
# Stage 2: runtime image (lean, split layers)
# =========================================
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    VENV_PATH=/var/www/FlaskApp/FlaskApp/venv

# Make apt more resilient + smaller indexes
RUN set -eux; \
  echo 'Acquire::Retries "5";' > /etc/apt/apt.conf.d/80-retries; \
  echo 'Acquire::http::No-Cache "true";' > /etc/apt/apt.conf.d/99no-cache; \
  echo 'Acquire::Languages "none";' > /etc/apt/apt.conf.d/99no-languages

# ---- Apache only (small layer)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends apache2; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /usr/share/doc/* /usr/share/man/* /usr/share/locale/*

# ---- mod_wsgi only (small layer)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends libapache2-mod-wsgi-py3; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /usr/share/doc/* /usr/share/man/* /usr/share/locale/*

# ---- Python runtime only (small layer)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends python3 python3-venv; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /usr/share/doc/* /usr/share/man/* /usr/share/locale/*

# ---- Core utilities (tiny layer)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends net-tools ca-certificates; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- Perl core (small)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- Perl libs batch 1 (small)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends libcgi-pm-perl libdata-stag-perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- Perl libs batch 2 (small)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends libfile-which-perl libio-string-perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- Perl libs batch 3 (small)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends libipc-run-perl libgo-perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- libgd-perl alone (can be chunky on 24.04)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends libgd-perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- Graphviz alone (largest piece, isolated)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends graphviz; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- Perl Graphviz bindings alone (separate from Graphviz)
RUN set -eux; \
  apt-get update || (sleep 3; apt-get update); \
  apt-get install -y --no-install-recommends libgraphviz-perl; \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# ---- (Optional) AWS CLI via pip in venv to keep layer isolated (off by default)
ARG INSTALL_AWSCLI=0
RUN if [ "$INSTALL_AWSCLI" = "1" ]; then \
      set -eux; \
      python3 -m venv /opt/awscli; \
      /opt/awscli/bin/pip install --upgrade pip setuptools wheel; \
      /opt/awscli/bin/pip install --no-cache-dir awscli; \
      ln -s /opt/awscli/bin/aws /usr/local/bin/aws; \
    fi

# ---- Bring in ONLY the built GO-TermFinder tree
COPY --from=builder /opt/gtf /opt/gtf

# Expose GO-TermFinder binaries/modules
ENV PATH="/opt/gtf/bin:${PATH}"
ENV PERL5LIB="/opt/gtf/lib/perl5:${PERL5LIB:-}"

# Copy web assets and Apache site config from your repo
COPY --from=builder /tmp/GoToolsDocker/www /var/www
COPY --from=builder /tmp/GoToolsDocker/FlaskApp.conf /etc/apache2/sites-available/FlaskApp.conf

# Ensure expected dirs & permissions
RUN set -eux; \
  install -d -m 1777 /var/www/tmp; \
  install -d -m 1777 /var/www/html; \
  install -d         /var/www/FlaskApp/FlaskApp; \
  install -d         /var/www/FlaskApp/FlaskApp/static

# Enable site
RUN set -eux; \
  a2enmod wsgi; \
  a2ensite FlaskApp; \
  a2dissite 000-default

# Create venv & install just what you need
RUN set -eux; \
  python3 -m venv "${VENV_PATH}"; \
  "${VENV_PATH}/bin/pip" install --upgrade pip setuptools wheel; \
  "${VENV_PATH}/bin/pip" install --no-cache-dir Flask flask-cors boto3; \
  rm -rf /root/.cache/pip

# Small healthcheck
HEALTHCHECK --interval=30s --timeout=5s --retries=5 CMD curl -fsS http://localhost/ || exit 1

EXPOSE 80
CMD ["apachectl", "-D", "FOREGROUND"]
