FROM ubuntu:latest

COPY . /server

WORKDIR /server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip python3-venv build-essential python3-dev && \
    python3 -m venv venv && \
    . venv/bin/activate && \
    pip install umbridge cuqipy

CMD . venv/bin/activate && python3 server.py
