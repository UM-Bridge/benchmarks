FROM mpioperator/openmpi-builder
#TODO: Separate builder and runtime images

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y build-essential git cmake pkg-config libsuitesparse-dev libarpack2-dev libarpack++2-dev openmpi-bin

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
    git clone https://gitlab.dune-project.org/anne.reinarz/dune-composites.git

COPY release.opts release.opts

RUN ./dune-common/bin/dunecontrol --opts=release.opts --builddir=$PWD/build all

WORKDIR /dune/build/dune-composites/src

RUN apt install -y python3-pip && pip3 install umbridge

COPY server.py server.py

CMD python3 server.py
