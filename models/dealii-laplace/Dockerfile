FROM dealii/dealii

COPY . /laplace/

RUN cd / && sudo git clone https://github.com/UM-Bridge/umbridge.git

RUN sudo mkdir /laplace/build &&\
    cd /laplace/build &&\
    sudo cmake -DCMAKE_EXE_LINKER_FLAGS="-lpthread -lrt" -DCMAKE_CXX_FLAGS=-I\ /umbridge/lib .. &&\
    sudo make

CMD /laplace/build/mcmc-laplace
