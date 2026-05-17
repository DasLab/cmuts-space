FROM python:3.11-slim-bookworm

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    git \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libhdf5-dev \
    libhts-dev \
    libomp-dev \
    bowtie2 \
    samtools \
    fontconfig \
    fonts-urw-base35 \
    && rm -rf /var/lib/apt/lists/*

# Debian bookworm ships cmake 3.25, cmuts needs 3.29+
RUN pip install --no-cache-dir cmake

# Clone and build cmuts at a pinned commit.
# This SHA is automatically updated by the "Sync HF Space" workflow in
# github.com/hmblair/cmuts on every push to master — do not pin manually.
ARG CMUTS_SHA=f97a05dc6672b1122ecf80be15a4221ff67627aa
RUN git clone --recurse-submodules https://github.com/hmblair/cmuts.git /cmuts && \
    cd /cmuts && git checkout $CMUTS_SHA
WORKDIR /cmuts
RUN ./configure

# Make htscodecs libs available
ENV PATH="/cmuts/bin:$PATH"
RUN cp /cmuts/htscodecs/lib/*.so* /usr/local/lib/ && ldconfig

# Patch Helvetica -> Nimbus Sans (available from fonts-urw-base35)
RUN sed -i 's/font.family.*=.*"Helvetica"/font.family"] = "Nimbus Sans"/' \
    /cmuts/src/python/cmuts/visualize/plotting.py && \
    sed -i 's/FONT_FAMILY = "Helvetica"/FONT_FAMILY = "Nimbus Sans"/' \
    /cmuts/src/python/cmuts/visualize/plotly.py && \
    fc-cache -f && \
    python3 -c "import matplotlib.font_manager; matplotlib.font_manager._load_fontmanager(try_read_cache=False)"

# Install Python web deps
RUN pip install --no-cache-dir \
    "fastapi>=0.110" \
    "uvicorn[standard]>=0.27" \
    "jinja2>=3.1" \
    "python-multipart>=0.0.9" \
    plotly "kaleido==0.2.1" h5py

# Clean up build artifacts
RUN rm -rf /cmuts/build

# Copy the app, templates, static assets, and example data
COPY app.py pipeline.py /app/
COPY templates /app/templates
COPY static /app/static
COPY examples /app/examples
WORKDIR /app

EXPOSE 7860

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH="/cmuts/bin:/home/user/.local/bin:$PATH"

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7860", "--workers", "1"]
