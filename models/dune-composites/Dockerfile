FROM mpioperator/openmpi
#TODO: Separate builder and runtime images

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y build-essential git cmake pkg-config libsuitesparse-dev libarpack2-dev libparpack2-dev libarpack++2-dev libopenmpi-dev g++

# Create directory and change working directory
RUN mkdir -p /dune
WORKDIR /dune

RUN git clone https://gitlab.dune-project.org/core/dune-common.git && \
    git clone https://gitlab.dune-project.org/core/dune-geometry.git && \
    git clone https://gitlab.dune-project.org/staging/dune-uggrid.git && \
    git clone https://gitlab.dune-project.org/core/dune-grid.git && \
    git clone https://gitlab.dune-project.org/core/dune-localfunctions.git && \
    git clone https://gitlab.dune-project.org/staging/dune-typetree.git && \
    git clone https://gitlab.dune-project.org/core/dune-istl.git && \
    git clone https://gitlab.dune-project.org/staging/dune-functions.git && \
    git clone https://gitlab.dune-project.org/pdelab/dune-pdelab.git && \
    git clone --branch=gcpscaling https://gitlab.dune-project.org/anne.reinarz/dune-composites.git

COPY release.opts release.opts

RUN ./dune-common/bin/dunecontrol --opts=release.opts --builddir=$PWD/build all

WORKDIR /dune/build/dune-composites/src

ENV OMP_NUM_THREADS=1

RUN apt install -y python3-pip python3-venv && \
    python3 -m venv venv && \
    . venv/bin/activate && \
    pip install umbridge

COPY server.py server.py

CMD . venv/bin/activate && python3 server.py
