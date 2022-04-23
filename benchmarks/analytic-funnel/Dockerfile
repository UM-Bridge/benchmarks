FROM ubuntu:latest

COPY . /server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-aiohttp python3-requests python3-numpy

CMD python3 /server/server.py
