FROM ubuntu:mantic

# Install fenics and python packages
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y bash fenics python3-pip
RUN pip3 install --break-system-packages umbridge

# Expose port for UM-BRIDGE server
EXPOSE 4242

RUN apt update

# Set up internal folders and files
WORKDIR /
COPY umbridge-server.py /
COPY cookiepde.py /

# Start Python UM-BRDIGE server
CMD python3 umbridge-server.py