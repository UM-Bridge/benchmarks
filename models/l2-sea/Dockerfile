FROM ubuntu:latest

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip python3-venv gfortran liblapack-dev libblas-dev libomp-dev git && \
    python3 -m venv venv && . venv/bin/activate && \
    pip install umbridge f90nml

RUN git clone https://github.com/MAORG-CNR-INM/NATO-AVT-331-L2-Sea-Benchmark.git

RUN cd NATO-AVT-331-L2-Sea-Benchmark && \
    make

# Symlink model output directory to /output for easier bind mounts
RUN mkdir -p /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000 && ln -s /NATO-AVT-331-L2-Sea-Benchmark/examples/DTMB-5415/CPU000 /output

COPY umbridge-server.py /

CMD ulimit -s unlimited && \
    . venv/bin/activate && \
    python3 umbridge-server.py
