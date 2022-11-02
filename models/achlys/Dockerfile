FROM achlys-ubuntu

RUN apt update && apt install -y inotify-tools && \
    rm -rf /var/lib/apt/lists/*

RUN apt install -y python3-pip && pip3 install umbridge

COPY server.py /server.py

COPY modify_input_file /opt/modify_input_file

CMD python3 /server.py
