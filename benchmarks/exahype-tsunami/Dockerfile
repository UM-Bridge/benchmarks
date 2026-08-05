FROM linusseelinger/model-exahype-tsunami

COPY . /benchmark-server

RUN apt update && \
    DEBIAN_FRONTEND="noninteractive" apt install -y build-essential zlib1g-dev \
    libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libsqlite3-dev \
    libreadline-dev libffi-dev curl libbz2-dev liblzma-dev

RUN curl -O https://www.python.org/ftp/python/3.8.2/Python-3.8.2.tar.xz

RUN tar -xf Python-3.8.2.tar.xz

RUN cd Python-3.8.2 && \
    ./configure --enable-optimizations --enable-loadable-sqlite-extensions && \
    make && make altinstall

RUN pip3.8 install --upgrade pip

RUN pip3.8 install scipy umbridge

CMD /server/server & python3.8 /benchmark-server/benchmark-server.py
