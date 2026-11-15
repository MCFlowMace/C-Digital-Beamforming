FROM nvidia/cuda:11.1.1-devel-ubuntu20.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential gcc-9 g++-9 \
        libarmadillo-dev \
        python3 python3-dev python-is-python3 \
        python3-setuptools python3-distutils cython3 \
        python3-numpy python3-matplotlib python3-scipy \
        python3-seaborn python3-statsmodels \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ENV CC=gcc-9 \
    CXX=g++-9 \
    CUDAHOSTCXX=g++-9

ENV MPLBACKEND=Agg

WORKDIR /project
CMD ["bash"]