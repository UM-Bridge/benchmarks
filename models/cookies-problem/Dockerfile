FROM ubuntu:mantic

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y autoconf automake autoconf-archive bash python3-pip pbzip2 wget clang-17 git cmake curl zip pkg-config gfortran lld libomp5-17 libomp-17-dev
RUN pip3 install --break-system-packages umbridge


RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list

RUN apt update
RUN apt install -y intel-basekit
# RUN . /opt/intel/oneapi/setvars.sh

RUN git clone https://gitlab.com/max.martinelli/igatools.git && \
  cd igatools && \
  git submodule init && \
  git submodule update && \
  git checkout fc6135ca63ed9062eefedfec5fda19eb3e3f4769 && \
  cd vcpkg && ./bootstrap-vcpkg.sh

#  git checkout 576316cfcec2be3212cba7a95f2b2df60a3f64c2 && \
#  git checkout a335daa3bf480cb6e57e2ea8522e3817e841a66d && \


WORKDIR /build_igatools
RUN . /opt/intel/oneapi/setvars.sh && cd /build_igatools && \
  cmake /igatools \
  -G"Unix Makefiles" \
  -DIGATOOLS_COMPONENT_DOCUMENTATION=OFF \
  -DIGATOOLS_WITH_GMP=ON \
  -DIGATOOLS_WITH_MKL=ON \
  -DMKL_DIR=/opt/intel/oneapi/mkl/latest \
  -DIGATOOLS_WITH_QUADRUPLE_PRECISION=OFF \
  -DIGATOOLS_WITH_SERIALIZATION=OFF \
  -DIGATOOLS_WITH_SUITESPARSE=OFF \
  -DIGATOOLS_WITH_SUPERLU=OFF \
  -DIGATOOLS_WITH_TRILINOS=OFF \
  -DTBB_DIR=/opt/intel/oneapi/tbb/latest \
  -DTBB_LIBRARY_DIR=/opt/intel/oneapi/tbb/latest/lib/intel64/gcc4.8/ \
  -DTRILINOS_DIR=/opt/trilinos/current \
  -DIGATOOLS_WITH_THREADS=ON \
  -Wno-dev \
  -DCMAKE_INSTALL_PREFIX=igatools \
  -DCMAKE_BUILD_TYPE=Release \
  -DIGATOOLS_STATIC_EXECUTABLE=OFF \
  -DCMAKE_C_COMPILER=clang-17 \
  -DCMAKE_CXX_COMPILER=clang++-17 \
  -DBUILD_SHARED_LIBS=ON

RUN make -j4 && \
  make setup_tests && \
  cd tests/ && \
  make poisson_lorenzo.release

WORKDIR /

COPY umbridge-server.py /

CMD python3 umbridge-server.py
