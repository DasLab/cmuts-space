FROM python:3.12-slim-bookworm

# System dependencies: the cmuts build chain plus the programs cmuts align calls.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    git \
    ca-certificates \
    libhts-dev \
    libhdf5-dev \
    zlib1g-dev \
    minimap2 \
    samtools \
    fastp \
    && rm -rf /var/lib/apt/lists/*

# Build and install cmuts at a pinned commit. The "Pin space to this commit"
# workflow in DasLab/cmuts updates this line on every push to its main branch,
# so do not pin it by hand.
ARG CMUTS_SHA=a1863a6086666b554840b1f27027692b425f973b
RUN git clone https://github.com/DasLab/cmuts.git /cmuts && \
    cd /cmuts && git checkout $CMUTS_SHA && \
    make && make install PREFIX=/usr/local && \
    rm -rf /cmuts

# Python dependencies for the app and for cmuts plot, read from the one
# list pyproject.toml holds.
COPY pyproject.toml /app/pyproject.toml
RUN pip install --no-cache-dir uv && \
    uv pip install --system --no-cache -r /app/pyproject.toml

COPY app.py pipeline.py options.py report.py job_description.py /app/
COPY templates /app/templates
COPY static /app/static
COPY examples /app/examples
WORKDIR /app

EXPOSE 7860

RUN useradd -m -u 1000 user
USER user

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7860", "--workers", "1"]
