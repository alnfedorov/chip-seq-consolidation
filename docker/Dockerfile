FROM ubuntu:bionic

ARG PYTHON=python3.8

# General tools
ARG FASTQC_VERSION=0.11.9
ARG CHIPS_VERSION=2.0
ARG HTSLIB_VERSION=1.8
ARG BOWTIE2_VERSION=2.4.1
ARG SAMBAMBA_VERSION=0.7.1
ARG PICARD_VERSION=2.22.1
ARG BWA_VERSION=0.7.17
ARG BEDTOOLS_VERSION=2.29.2

# Peak calling tools
ARG MACS2_VERSION=2.2.6
ARG EPIC2_VERSION=0.0.41
ARG PEAKRANGER_VERSION=1.18
ARG PEPR_VERSION=1.1.24
ARG MUSIC_COMMIT=b1caecdb164b1ab80acccb9463abe2526a56f69f
ARG IDR_VERSION=2.0.3

# dependencies and tools
RUN apt-get update && apt-get install -y \
	# general
	wget curl unzip ${PYTHON} ${PYTHON}-dev \
	# fastqc
	default-jre perl \
	# multiqc
	python3-pip \
	# chips
	autoconf-archive build-essential pkg-config libtool autoconf automake zlib1g-dev \
	# htslib
	libbz2-dev liblzma-dev \
	# picard
	openjdk-8-jdk \
	# samtools
	samtools \
	# matplotlib
	libfreetype6-dev

# general
RUN ${PYTHON} -m pip install Cython jupyter

RUN ${PYTHON} -m pip install \
	# tools
	multiqc \
	# peak calling deps
	pysam numpy \
    # used frameworks
    joblib pybedtools seaborn

WORKDIR /tools

# bedtools
RUN curl -SL https://github.com/arq5x/bedtools2/archive/v${BEDTOOLS_VERSION}.tar.gz | tar -xzv && \
	cd bedtools2-${BEDTOOLS_VERSION} && \
	make install && \
	cd .. && rm -rf bedtools2-${BEDTOOLS_VERSION}

# sambamba
RUN wget https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-static.gz && \
	gunzip sambamba-${SAMBAMBA_VERSION}-linux-static.gz && \
	rm -rf sambamba-${SAMBAMBA_VERSION}-linux-static.gz && \
	chmod +x sambamba-${SAMBAMBA_VERSION}-linux-static && \
	mv sambamba-${SAMBAMBA_VERSION}-linux-static /usr/bin/sambamba

# fastqc
# There is some problem with the server, wget stucks in the middle of the downloading.
RUN wget -c --tries=0 --read-timeout=20 http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && \
    unzip fastqc_v${FASTQC_VERSION}.zip && \
    rm fastqc_v${FASTQC_VERSION}.zip && \
    chmod +x FastQC/fastqc && \
    ln -s $(pwd)/FastQC/fastqc /usr/local/bin/fastqc

# chipmunk/tulip/chips - this tool has a lot of names
# depends on the htslib

RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
	tar -xjvf htslib-${HTSLIB_VERSION}.tar.bz2 && \
	rm htslib-${HTSLIB_VERSION}.tar.bz2 && \
	cd htslib-${HTSLIB_VERSION} && \
	./configure && \
	make && \
	make install

ENV PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/lib/pkgconfig
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/lib

RUN	wget https://codeload.github.com/gymreklab/chips/tar.gz/v${CHIPS_VERSION} && \
	tar -xzvf v${CHIPS_VERSION} && \
	rm v${CHIPS_VERSION} && \
	cd chips-${CHIPS_VERSION} && \
	# hard-patch sources
	sed -i '3i #include <random>' src/sequencer.h && \
	sed -i '3i #include <random>' src/pulldown.h && \
	autoreconf -if && ./configure && make install && \
	cd .. && rm -rf chips-${CHIPS_VERSION}

# picard
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar && \
	echo '#!/usr/bin/java -jar' > /usr/bin/picard && \
	cat picard.jar >> /usr/bin/picard && \
	rm -rf picard.jar && \
	chmod +x /usr/bin/picard

# bowtie2 aligner
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
	unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
	rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    ln -s $(pwd)/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2-build /usr/local/bin/bowtie2-build && \
    ln -s $(pwd)/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2 /usr/local/bin/bowtie2

# bwa aligner
RUN wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 && \
	tar -xjvf bwa-${BWA_VERSION}.tar.bz2 && \
	cd bwa-${BWA_VERSION} && make && cp bwa /usr/local/bin/bwa && cd .. && \
	rm -rf bwa-${BWA_VERSION}.tar.bz2 bwa-${BWA_VERSION}

# Peak-calling tools
RUN ${PYTHON} -m pip install macs2==${MACS2_VERSION} \
						     epic2==${EPIC2_VERSION} \
						   	 pepr==${PEPR_VERSION}

RUN wget https://excellmedia.dl.sourceforge.net/project/ranger/PeakRanger-${PEAKRANGER_VERSION}-Linux-x86_64.zip && \
	unzip -d PeakRanger PeakRanger-${PEAKRANGER_VERSION}-Linux-x86_64.zip && \
	rm -rf PeakRanger-${PEAKRANGER_VERSION}-Linux-x86_64.zip && \
	ln -s $(pwd)/PeakRanger/bin/peakranger /usr/local/bin/peakranger

RUN wget https://codeload.github.com/gersteinlab/MUSIC/zip/${MUSIC_COMMIT} && \
	unzip -d MUSIC ${MUSIC_COMMIT} && \
	rm -rf ${MUSIC_COMMIT} && \
	cd MUSIC/MUSIC-${MUSIC_COMMIT} && \
	make && \
	ln -s $(pwd)/bin/MUSIC /usr/local/bin/MUSIC && \
	chmod +x $(pwd)/bin/run_MUSIC.csh && \
	ln -s $(pwd)/bin/run_MUSIC.csh /usr/local/bin/run_MUSIC.csh

RUN wget https://github.com/nboley/idr/archive/${IDR_VERSION}.zip && \
	unzip ${IDR_VERSION}.zip && \
	cd idr-${IDR_VERSION} && ${PYTHON} setup.py install && cd .. && \
	rm -rf ${IDR_VERSION}.zip idr-${IDR_VERSION}

# SSH server
ARG public_key
EXPOSE 22

RUN apt-get update && apt-get install -y openssh-server \
    && mkdir -p /var/run/sshd \
    && mkdir -p ~/.ssh \
    && chmod 700 ~/.ssh \
    && sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin without-password/' /etc/ssh/sshd_config \
    && sed -i 's/#PasswordAuthentication yes/PasswordAuthentication no/' /etc/ssh/sshd_config \
    && echo ${public_key} >> ~/.ssh/authorized_keys \
    && echo "service ssh restart" >> ~/.bashrc  # hard-fix
