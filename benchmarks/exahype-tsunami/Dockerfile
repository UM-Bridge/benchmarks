FROM linusseelinger/model-exahype-tsunami

COPY . /benchmark-server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y python3-aiohttp python3-requests python3-scipy

CMD /server/server & python3 /benchmark-server/benchmark-server.py
