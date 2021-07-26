FROM python:3.6.9-slim-buster

ENV DEBIAN_FRONTEND noninteractive 
 
RUN apt-get update && apt-get install -y --no-install-recommends --no-install-suggests \ 
  git \
  vim \
  libqt4-network \
  libgomp1 \
  zlib1g \
  libstdc++6 \
  sudo \
  build-essential \
  libcurl4-openssl-dev

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

RUN echo "deb http://http.debian.net/debian buster main" > /etc/apt/sources.list.d/debian-unstable.list

RUN apt-get install -y --no-install-recommends \
    gfortran \
    pandoc

RUN pip install pip --upgrade
RUN pip install setuptools --upgrade
RUN pip install numpy mypy jupyter flask pytest
RUN pip install --force-reinstall 'pytest<=5.0.1'
RUN pip install -U git+https://github.com/bigbio/mzqc-pylib.git@v1.0.0
RUN pip install biopython pandas requests plotly-express toposort click
RUN pip install pyopenms==2.5.* 

#RUN pip install -e .  # devcontainer.json: "postCreateCommand": "pip install -e ."

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME 
RUN echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME
RUN echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME
    
#USER 1000  # "runArgs": ["-u", "1000"]
#extrawurst fuer vscode
#RUN mkdir -p /root