FROM ubuntu:latest

COPY umbridge-server.py /

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip gfortran liblapack-dev libblas-dev libomp-dev git && \
    pip install umbridge

RUN git clone https://github.com/MAORG-CNR-INM/NATO-AVT-331-L2-Sea-Benchmark.git

RUN cd NATO-AVT-331-L2-Sea-Benchmark && \
    make

CMD ulimit -s unlimited && \
    python3 umbridge-server.py
