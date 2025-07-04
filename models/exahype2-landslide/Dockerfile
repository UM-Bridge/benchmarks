FROM peanoframework/exahype2:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install umbridge

RUN cd /exahype2/applications/exahype2/shallow-water/tafjord-landslide/

COPY ./server.py /

RUN python3 tafjord-landslide.py --friction 100 -ns 0

CMD python3 /server.py
