FROM linusseelinger/model-achlys

COPY . /benchmark-server

CMD . venv/bin/activate && python3 /server.py & \
    . venv/bin/activate && python3 /benchmark-server/benchmark-server.py
