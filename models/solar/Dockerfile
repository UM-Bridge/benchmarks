FROM ubuntu:latest

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip git && \
    pip install umbridge

RUN git clone https://github.com/bbopt/solar.git

RUN cd solar/src && make

COPY server.py  /

CMD python3 server.py
