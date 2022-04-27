FROM mparno/muq:latest

COPY . /server

USER root
RUN apt update && \
    apt install -y python3-aiohttp python3-requests python3-numpy python3-h5py

USER muq-user

WORKDIR /server 
CMD python3 BeamModel.py
