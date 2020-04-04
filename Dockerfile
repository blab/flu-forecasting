# Start from a slim Debian image.
FROM debian:buster-slim
LABEL MAINTAINER="John Huddleston <huddlej@gmail.com>"

ENV CONDA_VERSION 4.7.12.1
ENV CONDA_MD5 81c773ff87af5cfac79ab862942ab6b3

# Tell Python not to recreate the bytecode files. Since this is a docker image,
# these will be recreated every time, writing them just uses unnecessary disk
# space.
ENV PYTHONDONTWRITEBYTECODE=true

# Setup miniconda after the official Dockerfile's approach:
# https://hub.docker.com/r/continuumio/miniconda3/dockerfile
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Use bash as the default shell.
SHELL ["/bin/bash", "-c"]

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
      wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-$CONDA_VERSION-Linux-x86_64.sh \
    && echo "${CONDA_MD5}  Miniconda3-$CONDA_VERSION-Linux-x86_64.sh" > miniconda.md5 \
    && if [ $(md5sum -c miniconda.md5 | awk '{print $2}') != "OK" ] ; then exit 1; fi \
    && mv Miniconda3-$CONDA_VERSION-Linux-x86_64.sh miniconda.sh \
    && /bin/bash ./miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh miniconda.md5 \
    && /opt/conda/bin/conda clean -tipsy \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate base" >> ~/.bashrc \
    && /opt/conda/bin/conda install --freeze-installed tini -y \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && /opt/conda/bin/conda clean -afy

# Create Nextstrain environment inside the base environment.
COPY envs/anaconda.python3.yaml /tmp
WORKDIR /tmp
RUN conda env update -n base -f anaconda.python3.yaml \
 && conda clean --all --yes --force-pkgs-dirs \
 && find /opt/conda/ -follow -type f -name '*.pyc' -delete

# Create a top-level directory for all Nextstrain-related tools.
WORKDIR /nextstrain

# Download the fauna repository.
RUN git clone https://github.com/nextstrain/fauna.git

# Copy files to run the forecasting pipeline.
WORKDIR /nextstrain/flu-forecasting
# COPY Snakefile .
# COPY config config/
# COPY data data/
# COPY dist dist/
# COPY envs envs/
# COPY rules rules/
# COPY scripts scripts/
# COPY src src/

ENTRYPOINT ["/opt/conda/bin/tini", "-g", "--", "/bin/bash"]
