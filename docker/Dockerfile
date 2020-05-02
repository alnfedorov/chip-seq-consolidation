FROM nvidia/cuda:10.2-cudnn7-devel

ARG PYTHON=python3.8

# General tools
ARG BEDTOOLS_VERSION=2.29.2

# dependencies and tools
RUN apt-get update && apt-get install -y \
	# general
	wget curl unzip python ${PYTHON} ${PYTHON}-dev \
	# multiqc
	python3-pip \
	# htslib
	libbz2-dev liblzma-dev \
	# matplotlib
	libfreetype6-dev

# general
RUN ${PYTHON} -m pip install Cython numpy jupyter

RUN ${PYTHON} -m pip install \
	# tools
	multiqc \
	# peak calling deps
	pysam \
    # used frameworks
    joblib pybedtools matplotlib seaborn pyBigWig torch pytorch-ignite tqdm numba scikit-learn

WORKDIR /tools

# bedtools
RUN curl -SL https://github.com/arq5x/bedtools2/archive/v${BEDTOOLS_VERSION}.tar.gz | tar -xzv && \
	cd bedtools2-${BEDTOOLS_VERSION} && \
	make install && \
	cd .. && rm -rf bedtools2-${BEDTOOLS_VERSION}

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
