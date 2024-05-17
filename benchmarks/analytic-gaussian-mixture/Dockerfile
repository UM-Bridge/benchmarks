FROM ubuntu:latest

COPY . /server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip python3-venv && \
    python3 -m venv venv && \
    . venv/bin/activate && \
    pip install umbridge numpy scipy

CMD . venv/bin/activate && python3 /server/server.py
