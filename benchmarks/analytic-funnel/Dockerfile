FROM ubuntu:latest

COPY . /server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip

RUN pip3 install umbridge numpy

CMD python3 /server/server.py
