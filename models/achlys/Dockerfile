FROM achlys-ubuntu

RUN apt update && apt install -y inotify-tools && \
    rm -rf /var/lib/apt/lists/*

RUN apt install -y python3-pip && pip3 install umbridge

RUN cd /opt && git clone --branch=umbridge https://github.com/digiLab-ai/achlys-uq.git

COPY server.py /server.py

CMD python3 /server.py
