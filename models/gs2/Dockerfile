FROM ubuntu:20.04
LABEL maintainer="william.hornsby@ukaea.uk"
LABEL version="1.0"
LABEL description="This is Docker Image for running GS2 with GNU compilers on \
Ubuntu 20.04. This uses MPICH instead of OpenMPI to avoid issues with running \
as root."

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
	python3.8-venv \
	python3-pip \
	mpich \
	libnetcdf-dev \
	libnetcdff-dev \
	libfftw3-dev \
        python-is-python3

RUN echo 'alias python=python3' >> $HOME/.bashrc
RUN echo 'alias pip=pip3' >> $HOME/.bashrc

RUN echo 'export GK_SYSTEM=gnu_ubuntu' >> $HOME/.bashrc

RUN pip install ford netCDF4==1.6.2

RUN git clone https://whornsby@bitbucket.org/gyrokinetics/gs2.git /usr/gs2

WORKDIR /usr/gs2

RUN cd /usr/gs2

ENV GK_SYSTEM=gnu_ubuntu

RUN echo ${GK_SYSTEM}

RUN make -IMakefiles depend

RUN make -IMakefiles gs2

#RUN make -IMakefiles tests
