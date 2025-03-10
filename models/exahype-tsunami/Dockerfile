FROM chun9l/openmpi-builder:0.3.0 as builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y python3 git cmake wget curl zlib1g pkg-config libhdf5-openmpi-dev libc6-dev libnuma-dev doxygen libyaml-cpp-dev gfortran libtbb-dev libstdc++-8-dev liblua5.3-dev

RUN cd / && mkdir dependencies && \
    cd /dependencies/ && git clone https://github.com/uphoffc/ImpalaJIT.git && \
    cd ImpalaJIT/ && \
    mkdir build && cd build &&\
    CC=mpicc CXX=mpicxx cmake -DSHARED_LIB=1 .. && make -j4 && make install

RUN cd /dependencies && wget https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/netcdf-c-4.8.1.tar.gz &&\
    tar -xvzf netcdf-c-4.8.1.tar.gz && cd netcdf-c-4.8.1 &&\
    CC=mpicc CXX=mpicxx CPPFLAGS=-I/usr/include/hdf5/openmpi/ LDFLAGS=-L/usr/lib/$(arch)-linux-gnu/hdf5/openmpi/ ./configure --enable-parallel &&\
    #CC=mpicc CXX=mpicxx ./configure -prefix=/dependencies/netcdf-c-4.8.1/lib --enable-parallel --disable-dap  &&\
    make install

ENV PKG_CONFIG_PATH="/dependencies/ImpalaJIT/build/:/dependencies/netcf-c-4.8.1/lib/lib/pkgconfig:/usr/lib/$(arch)-linux-gnu/pkgconfig"
ENV CMAKE_PREFIX_PATH=/dependencies/netcdf-c-4.8.1/lib/

RUN cd /dependencies && git clone --recursive https://github.com/TUM-I5/ASAGI.git && \
    cd ASAGI/ && \
    mkdir build && cd build && \
    NetCDF_LIBRARY=/dependencies/netcdf-c-4.8.1/lib cmake -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
               -DTHREADSAFE=ON -DNONUMA=ON -DSHARED_LIB=1 -DCMAKE_CXX_FLAGS="-I/dependencies/netcdf-c-4.8.1/include/" ..  && \
    make -j4 && make install

ENV COMPILER_LFLAGS=" -L/dependencies/ImpalaJIT/build/ -limpalajit -L/dependencies/netcdf-c-4.8.1/lib -lnetcdf  -L/usr/lib/$(arch)-linux-gnu/ -L/usr/local/lib/ -lyaml-cpp -lasagi -L/dependencies/easi/build -leasi -lpthread"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib/:/dependencies/ImpalaJIT/build"

RUN cd /dependencies/ && git clone https://github.com/SeisSol/easi.git && \
    cd easi && git checkout 647ff0863a83c011e56c8eb6b4fc522fa355c153 &&\
    mkdir build && cd build && \
    CC=mpicc CXX=mpicxx cmake -DPKG_CONFIG_USE_CMAKE_PREFIX_PATH=TRUE -DCMAKE_PREFIX_PATH=/dependencies/ -DCMAKE_INSTALL_PREFIX=/dependencies/easi/build/ .. && \
    make -j4

RUN cd / && git   clone https://github.com/annereinarz/ExaHyPE-Tsunami.git
WORKDIR /ExaHyPE-Tsunami

ENV COMPILER_LFLAGS="-DMPI2 -lmpi_cxx -lnuma -fopenmp -lstdc++fs -llua5.3 -L/dependencies/ImpalaJIT/build/ -limpalajit -L/dependencies/netcdf-c-4.8.1/lib/lib -lnetcdf -L/usr/lib/$(arch)-linux-gnu/ -lyaml-cpp -L/usr/local/lib -lasagi -L/dependencies/easi/build -leasi -lpthread -ltbb"
ENV COMPILER_CFLAGS=" -I/dependencies/easi/include/"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib:/dependencies/easi/build"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/dependencies/netcdf-c-4.8.1/lib/lib/"
ENV SHAREDMEM=None
ENV COMPILER=GNU
#ENV TBB_INC=/usr/include/tbb
#ENV TBB_SHLIB=/usr/lib/$(arch)-linux-gnu/libtbb.so
ENV EXAHYPE_CC=mpicc
ENV EXAHYPE_FC=mpicc++
ENV PROJECT_CFLAGS="-DUSE_ASAGI"

RUN cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/ && \
    ../../Toolkit/toolkit.sh SWE_asagi_limited_l0.exahype2 && \
    ../../Toolkit/toolkit.sh SWE_asagi_limited_l1.exahype2 && \
    ../../Toolkit/toolkit.sh SWE_asagi_limited_l2.exahype2 && \
    cd SWE_asagi_limited_l0    && make -j4 && mkdir Probes  &&\
    cd ../SWE_asagi_limited_l1 && make -j4 && mkdir Probes  && \
    cd ../SWE_asagi_limited_l2 && make -j4 && mkdir Probes

# Default verbosity level... need cleaner solution
RUN cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l0 && cp exahype_debug.log-filter exahype.log-filter && \
    cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l1 && cp exahype_debug.log-filter exahype.log-filter && \
    cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l2 && cp exahype_debug.log-filter exahype.log-filter

RUN cd / && git clone https://github.com/UM-Bridge/umbridge.git

# Build model server
COPY . /server
RUN cd /server && \
    g++ -o server server.cpp -I /umbridge/lib -lpthread

FROM mpioperator/openmpi:0.3.0

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y python3 libyaml-cpp0.6 libtbb2 libhdf5-openmpi-103 liblua5.3-dev && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /server/server /server/server
COPY --from=builder /ExaHyPE-Tsunami /ExaHyPE-Tsunami
COPY --from=builder /dependencies /dependencies
COPY --from=builder /usr/local /usr/local

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib/:/dependencies/ImpalaJIT/build"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib:/dependencies/easi/build"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/dependencies/netcdf-c-4.8.1/lib/lib/"

#RUN cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l0 && mkdir Probes && mkdir vtk-output && \
#    cd ../SWE_asagi_limited_l1 && mkdir Probes && mkdir vtk-output && \
#    cd ../SWE_asagi_limited_l2 && mkdir Probes && mkdir vtk-output

WORKDIR /ExaHyPE-Tsunami

ENV PORT=4242
ENV RANKS=1
ENV SHARED_DIR=/shared/comm
RUN mkdir /output /shared
CMD /server/server
