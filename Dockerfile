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

# Clone and build cmuts at a pinned commit (updated by CI)
ARG CMUTS_SHA=cfc447be42ddd42db87452e4836c88ce7aa170a4
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
    fc-cache -f && \
    python3 -c "import matplotlib.font_manager; matplotlib.font_manager._load_fontmanager(try_read_cache=False)"

# Install Gradio
RUN pip install --no-cache-dir gradio plotly h5py uvicorn

# Clean up build artifacts
RUN rm -rf /cmuts/build

# Copy the app and example data
COPY app.py /app/app.py
COPY examples /app/examples
WORKDIR /app

EXPOSE 7860

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH="/cmuts/bin:/home/user/.local/bin:$PATH"

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7860"]
