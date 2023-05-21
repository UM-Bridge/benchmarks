FROM ubuntu:latest

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip

RUN pip3 install umbridge cuqipy

COPY . /server

WORKDIR /server

CMD python3 server.py
