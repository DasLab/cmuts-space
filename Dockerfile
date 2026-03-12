FROM python:3.11-slim-bookworm

# System dependencies for cmuts-core, bowtie2, samtools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pkg-config \
    git \
    wget \
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
    && rm -rf /var/lib/apt/lists/*

# Clone cmuts with submodules
RUN git clone --recurse-submodules https://github.com/hmblair/cmuts.git /cmuts
WORKDIR /cmuts

# Build htscodecs (bundled submodule)
RUN cd htscodecs && \
    autoreconf -i && \
    ./configure --prefix=/cmuts/htscodecs && \
    make -j$(nproc) && \
    make install && \
    ldconfig

# Build cmuts-core
RUN mkdir -p build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc)

# Install binaries to PATH
RUN cp /cmuts/bin/cmuts-core /usr/local/bin/ && \
    cp /cmuts/bin/cmuts /usr/local/bin/ && \
    cp /cmuts/bin/cmuts-align /usr/local/bin/ && \
    cp /cmuts/htscodecs/lib/*.so* /usr/local/lib/ && \
    ldconfig

# Install Python package + Gradio
RUN pip install --no-cache-dir /cmuts gradio

# Clean up build artifacts
RUN rm -rf /cmuts/build

# Copy the app
COPY app.py /app/app.py
WORKDIR /app

# HF Spaces expects port 7860
EXPOSE 7860

# Run as non-root user (HF Spaces requirement)
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH="/home/user/.local/bin:$PATH"

CMD ["python", "app.py"]
