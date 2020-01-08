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

RUN pip install pip --upgrade
RUN pip install setuptools --upgrade
RUN pip install numpy mypy pyopenms==2.3.*
RUN pip install pronto pytest jupyter rpy2 biopython flask 
RUN pip install --force-reinstall 'pytest<=5.0.1'
RUN pip install -U git+https://github.com/bigbio/mzqc-pylib.git#egg=mzqc-pylib

COPY . /
RUN pip install -U -e .
