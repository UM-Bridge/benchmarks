FROM ubuntu:25.04

ENV DEBIAN_FRONTEND=noninteractive \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN apt-get update && apt-get install -y build-essential vim git-lfs openmpi-bin libopenmpi-dev cmake python3 python3-pip python3-venv libnetcdf-dev libxrender-dev

RUN git clone https://gitlab.lrz.de/hpcsoftware/Peano.git

WORKDIR /Peano

RUN git switch muc/exahype && mkdir build

RUN cd build && cmake ..

RUN cd build && make -j10

RUN python3 -m venv venv

RUN . venv/bin/activate && pip3 install -r requirements.txt && pip3 install -e . && pip3 install umbridge

RUN . venv/bin/activate && cd applications/exahype2/shallow-water/tohoku-tsunami && python3 tohoku-tsunami.py
