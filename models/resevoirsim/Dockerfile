FROM ubuntu:22.04

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y  build-essential git gfortran python3 python3-six python-six flex bison python3-pip python3-venv && \
	python3 -m venv venv && . venv/bin/activate && \
	pip3 install umbridge numpy scipy pandas && \
	rm -rf /var/lib/apt/lists/*

RUN git clone -b release https://bitbucket.org/opengosim/pflotran_ogs_1.8.git

RUN mkdir petsc && git clone https://gitlab.com/petsc/petsc.git petsc 

RUN cd petsc && git checkout v3.19.1

ENV PETSC_DIR=/petsc
ENV PETSC_ARCH=ubuntu-opt

RUN cd petsc && ./configure --download-mpich=yes --download-hdf5=yes --with-hdf5-fortran-bindings=yes --download-fblaslapack=yes --download-ptscotch=yes --download-hypre=yes --with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 && make PETSC_DIR=/petsc PETSC_ARCH=ubuntu-opt all && make PETSC_DIR=/petsc PETSC_ARCH=ubuntu-opt check

ENV USER=root

RUN cd pflotran_ogs_1.8/src/pflotran && make -j4 pflotran modern_gfortran=1

COPY server.py /server.py
COPY coarse_model /coarse_model

RUN mkdir /output

ENV PFT_PATH=/
RUN sed -i 's|^PFT_PATH=.*|PFT_PATH=/pflotran_ogs_1.8|' /pflotran_ogs_1.8/macros/pft1.8
ENV PORT=4243

CMD . venv/bin/activate && cd coarse_model && python3 /server.py
