FROM mpioperator/openmpi
LABEL maintainer="william.hornsby@ukaea.uk"
LABEL version="1.0"
LABEL description="This is Docker Image for running GS2 with GNU compilers on \
Debian. This uses OpenMPI."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \ 
        && apt install -y \
        gfortran \
        g++ \
        make \
        wget \
        git \
        vim \
        cmake \
        rsync \
        curl \
        python3 \
        python3-venv \
        python3-pip \
        libnetcdf-dev \
        libnetcdff-dev \
        netcdf-bin \
        libfftw3-dev \
        libfftw3-mpi-dev \
        python-is-python3 \
        openmpi-bin \
        openmpi-doc \
        libopenmpi-dev

RUN echo 'alias python=python3' >> $HOME/.bashrc
RUN echo 'alias pip=pip3' >> $HOME/.bashrc

RUN echo 'export GK_SYSTEM=gnu_ubuntu' >> $HOME/.bashrc

RUN git clone --recurse-submodules https://bitbucket.org/gyrokinetics/gs2.git /usr/gs2

WORKDIR /usr/gs2

RUN cd /usr/gs2

ENV GK_SYSTEM=gnu_ubuntu

RUN echo ${GK_SYSTEM}

RUN make -IMakefiles gs2
#RUN make -IMakefiles tests

RUN python -m venv venv && . venv/bin/activate && \
    pip install umbridge pyrokinetics 

COPY server.py fast.in .

CMD . venv/bin/activate && python server.py
