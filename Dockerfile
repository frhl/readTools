# Use an official base image, e.g., Ubuntu
FROM ubuntu:20.04

# Set maintainer label
LABEL maintainer="flassen@well.ox.ac.uk"

# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhts-dev \
    git

# Clone htslib
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git

# Install htslib
WORKDIR htslib
RUN make
RUN make install
WORKDIR ..

# copy scripts
WORKDIR app
COPY makefile makefile
COPY candiates.cpp candidates.cpp
COPY phaseReads.cpp phaseReads.cpp
RUN make

# move to folder in PATH
RUN mv candidates /usr/local/bin/.
RUN mv phaseReads /usr/local/bin/.

# Set default command to R when the container starts
#CMD ["bash"]


