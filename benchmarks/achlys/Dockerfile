FROM linusseelinger/model-achlys

COPY . /benchmark-server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip python3-venv

RUN python3 -m venv venv && . venv/bin/activate && \
    pip3 install umbridge scipy

CMD . venv/bin/activate && \
    /server/server & python3 /benchmark-server/benchmark-server.py
