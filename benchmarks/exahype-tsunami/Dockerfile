FROM linusseelinger/model-exahype-tsunami

COPY . /benchmark-server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-pip python3-scipy

RUN pip3 install umbridge

CMD /server/server & python3 /benchmark-server/benchmark-server.py
