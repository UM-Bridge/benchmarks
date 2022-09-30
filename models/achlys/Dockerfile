FROM achlys-ubuntu

RUN apt update && apt install -y inotify-tools && \
    rm -rf /var/lib/apt/lists/*

RUN apt install -y python3-pip && pip3 install umbridge

RUN cd /opt && git clone https://github.com/digiLab-ai/achlys-uq.git

#RUN source achlys-uq/scripts/bashrc && \
#    ./achlys-uq/achlys-uq/run_desorp

COPY server.py /server.py

RUN python3 /server.py
