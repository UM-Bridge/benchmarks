FROM ubuntu:latest

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y build-essential libssl-dev git cmake 

RUN cd / && git clone https://github.com/UM-Bridge/umbridge.git

RUN cd / && mkdir Code && cd Code && mkdir ocean && cd /

COPY core /Code/core 
COPY ocean/CMakeLists.txt  /Code/ocean/CMakeLists.txt
COPY ocean/Sources  /Code/ocean/Sources
COPY ocean/Data  /Code/ocean/Data

RUN cd /Code/ocean && \
    mkdir Build && cmake -S /Code/ocean -B Build && \ 
    cd /Code/ocean/Build/Sources && make 
    #&& ./WaveInOcean_um /Code/ocean/Data/ /Code/ocean/Build/Data/

#CMD ./server/minimal-server
CMD ./Code/ocean/Build/Sources/WaveInOcean_um /Code/ocean/Data/ /Code/ocean/Build/Data/