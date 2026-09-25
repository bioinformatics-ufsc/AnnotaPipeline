FROM ubuntu:20.04

LABEL maintainer="Laboratório de Bioinformática - UFSC"
LABEL description="AnnotaPipeline: Multi-Omics Protein Annotation Pipeline"

ENV DEBIAN_FRONTEND="noninteractive"
ENV TZ="Etc/UTC"

# Software versions
ENV IPR="5"
ENV HMMER_VERSION="3.4"
ENV IPRSCAN_VERSION="5.68-100.0"
ENV BLAST_VERSION="2.16.0"
ENV PERCOLATOR_VERSION="3-07-01"
ENV COMET_VERSION="2021010"
ENV XERCES_VERSION="3.2.4"

# Set timezone and install base packages & dependencies
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        build-essential \
        wget \
        curl \
        git \
        cmake \
        autoconf \
        zip \
        unzip \
        ca-certificates \
        rpm \
        python3 \
        python3-pip \
        python3-dev \
        python3-distutils \
        python3-biopython \
        openjdk-11-jre-headless \
        perl \
        perl-doc \
        libfile-which-perl \
        libparallel-forkmanager-perl \
        libyaml-perl \
        libdbd-mysql-perl \
        libpcre2-dev \
        libgomp1 \
        libgsl-dev \
        libboost-all-dev \
        libsuitesparse-dev \
        liblpsolve55-dev \
        libsqlite3-dev \
        libmysql++-dev \
        libbamtools-dev \
        samtools \
        libhts-dev \
        cdbfasta \
        diamond-aligner \
        libtokyocabinet-dev \
        zlib1g-dev \
        libbz2-dev && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get install -y --no-install-recommends g++-11 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Symlink python
RUN ln -sf /usr/bin/python3 /usr/bin/python

# --------------------------------
# 1. Install BLAST+
# --------------------------------
RUN mkdir -p /opt/blast && \
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz -O /tmp/blast.tar.gz && \
    tar -zxf /tmp/blast.tar.gz -C /opt/blast --strip-components=1 && \
    rm -f /tmp/blast.tar.gz
ENV PATH="/opt/blast/bin:$PATH"

# --------------------------------
# 2. Install InterProScan (Core)
# --------------------------------
RUN mkdir -p /opt/interproscan && \
    curl -sSL "ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/${IPR}/${IPRSCAN_VERSION}/alt/interproscan-core-${IPRSCAN_VERSION}.tar.gz" -o /tmp/iprscan.tar.gz && \
    tar -pxzf /tmp/iprscan.tar.gz -C /opt/interproscan --strip-components=1 && \
    rm -f /tmp/iprscan.tar.gz
ENV PATH="/opt/interproscan/bin:/opt/interproscan:$PATH"

# --------------------------------
# 3. Install HMMER
# --------------------------------
RUN mkdir -p /opt/hmmer && \
    wget -q "http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz" -O /tmp/hmmer.tar.gz && \
    tar -zxf /tmp/hmmer.tar.gz -C /opt/hmmer --strip-components=1 && \
    rm -f /tmp/hmmer.tar.gz && \
    cd /opt/hmmer && \
    ./configure --prefix=/opt/hmmer && \
    make -j$(nproc) && \
    make install
ENV PATH="/opt/hmmer/bin:$PATH"

# --------------------------------
# 4. Install AUGUSTUS
# --------------------------------
RUN git clone --depth 1 https://github.com/Gaius-Augustus/Augustus.git /opt/Augustus && \
    cd /opt/Augustus && \
    make augustus -j$(nproc)
ENV PATH="/opt/Augustus/bin:/opt/Augustus/scripts:$PATH"
ENV AUGUSTUS_CONFIG_PATH="/opt/Augustus/config/"

# --------------------------------
# 5. Install Kallisto
# --------------------------------
RUN git clone --depth 1 https://github.com/pachterlab/kallisto.git /opt/kallisto-src && \
    mkdir -p /opt/kallisto-src/build && \
    cd /opt/kallisto-src/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/opt/kallisto && \
    make -j$(nproc) && \
    make install && \
    ln -sf /opt/kallisto/bin/kallisto /usr/local/bin/kallisto && \
    rm -rf /opt/kallisto-src

# --------------------------------
# 6. Install Comet MS/MS
# --------------------------------
RUN mkdir -p /opt/comet && \
    wget -q "https://sourceforge.net/projects/comet-ms/files/comet_${COMET_VERSION}.zip" -O /tmp/comet.zip && \
    unzip -q /tmp/comet.zip -d /opt/comet && \
    mv /opt/comet/comet.${COMET_VERSION}.debian.exe /opt/comet/comet && \
    chmod 755 /opt/comet/comet && \
    ln -sf /opt/comet/comet /usr/local/bin/comet && \
    rm -rf /tmp/comet.zip /opt/comet/*.exe

# --------------------------------
# 7. Install Xerces-C and Percolator
# --------------------------------
RUN mkdir -p /opt/percolator/xerces && \
    wget -q --no-check-certificate "https://archive.apache.org/dist/xerces/c/3/sources/xerces-c-${XERCES_VERSION}.tar.gz" -O /tmp/xerces.tar.gz && \
    tar -zxf /tmp/xerces.tar.gz -C /opt/percolator/xerces --strip-components=1 && \
    rm -f /tmp/xerces.tar.gz && \
    cd /opt/percolator/xerces && \
    ./configure --prefix=/opt/percolator/xerces --disable-netaccessor-curl --disable-transcoder-icu && \
    make -j$(nproc) && \
    make install

RUN wget -q "https://github.com/percolator/percolator/archive/refs/tags/rel-${PERCOLATOR_VERSION}.tar.gz" -O /tmp/percolator.tar.gz && \
    mkdir -p /opt/percolator/src && \
    tar -zxf /tmp/percolator.tar.gz -C /opt/percolator/src --strip-components=1 && \
    rm -f /tmp/percolator.tar.gz && \
    mkdir -p /opt/percolator/src/build && \
    cd /opt/percolator/src/build && \
    cmake .. && \
    make -j$(nproc) && \
    ln -sf /opt/percolator/src/build/src/percolator /usr/local/bin/percolator

# --------------------------------
# 8. Install AnnotaPipeline
# --------------------------------
WORKDIR /opt/AnnotaPipeline
COPY . /opt/AnnotaPipeline
RUN pip3 install --no-cache-dir .

WORKDIR /data
VOLUME /data

CMD ["AnnotaPipeline", "--help"]