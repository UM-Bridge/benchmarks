FROM krosenf/emod_base:v2.21.0_Ubuntu22
WORKDIR /experiment
RUN pip install numpy umbridge
# add the template files
COPY . /experiment
CMD python3 /experiment/model.py

