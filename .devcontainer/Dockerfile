FROM mcr.microsoft.com/devcontainers/cpp:ubuntu-22.04

# Update and install Tier-1 C++ build tools and profilers
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
        software-properties-common \
        build-essential \
        ninja-build \
        cmake \
        gdb \
        valgrind \
        linux-tools-common \
        linux-tools-generic \
        libbenchmark-dev \
        libgtest-dev \
    && rm -rf /var/lib/apt/lists/*

# Install modern compilers (GCC-13 and Clang-17)
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y \
    && apt-get update \
    && apt-get install -y gcc-13 g++-13 clang-17 clangd-17

# Set GCC-13 as the default compiler
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100 \
    && update-alternatives --install /usr/bin/clangd clangd /usr/bin/clangd-17 100