FROM quay.io/fenicsproject/dev

USER root
RUN pip3 install --upgrade pip

USER fenics

RUN pip3 install hippylib && \
    pip3 install jupyter && \
    pip3 install matplotlib && \
    pip3 install h5py && \
    pip3 install pyyaml && \
    pip3 install seaborn==0.10.0 && \
    pip3 install statsmodels && \
    pip3 install aiohttp && \
    pip3 install requests

WORKDIR /home/fenics/

# Install MUQ
USER root
RUN apt update && apt install -y libboost-all-dev libhdf5-dev && \
    git clone --depth 1 https://bitbucket.org/mituq/muq2.git && \
    cd muq2/; mkdir build; cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/home/fenics/lib/muq -DMUQ_USE_MPI=OFF -DMUQ_USE_PYTHON=ON ../ && \
    make -j2 install && \
    cd /home/fenics/ && \
    rm -rf muq2
USER fenics

# Install hippylib-muq interface
RUN git clone https://github.com/hippylib/hippylib2muq.git

# Set environmental variables
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/home/fenics/lib/muq/lib" \
    PYTHONPATH="/home/fenics/lib/muq/python:/home/fenics/hippylib2muq"

COPY ppoisson_box.py ppoisson_box.yaml nonlinearPPoissonProblem.py minimization.py /home/fenics/app/

USER root

WORKDIR /home/fenics/app/
CMD ["python3 /home/fenics/app/ppoisson_box.py"]
