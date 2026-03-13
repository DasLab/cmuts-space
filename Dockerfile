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
ARG CMUTS_SHA=52675908d4ccc3f372fca327f939d3a94169170d
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
RUN pip install --no-cache-dir gradio

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

CMD ["python", "app.py"]
