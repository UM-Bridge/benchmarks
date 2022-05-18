FROM dealii/dealii

RUN sudo apt update && sudo apt install -y libeigen3-dev

COPY . /laplace/

RUN sudo mkdir /laplace/build && cd /laplace/build && sudo cmake .. && sudo make
