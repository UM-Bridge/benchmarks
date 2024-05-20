FROM linusseelinger/achlys-ubuntu

RUN apt update && apt install -y inotify-tools python3-pip python3-venv && \
    python3 -m venv venv && . venv/bin/activate && \
    pip3 install umbridge numpy scipy pandas && \
    rm -rf /var/lib/apt/lists/*

COPY server.py /server.py

COPY modify_input_file /opt/modify_input_file

CMD . venv/bin/activate && python3 /server.py
