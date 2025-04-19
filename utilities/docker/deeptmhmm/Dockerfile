# This image contains only dependencies for DeepTMHMM. It does not include a copy of DeepTMHMM
FROM python:3.8-slim
ENV TZ=Europe/London
WORKDIR /opt/tmhmm
COPY . .
RUN apt-get update -y && apt-get upgrade -y && apt-get install -y libhdf5-dev procps
RUN python3 -m pip install wheel Cython==0.29.37 pkgconfig==1.5.5 \
  && python3 -m pip install https://download.pytorch.org/whl/cu92/torch-1.5.0%2Bcu92-cp38-cp38-linux_x86_64.whl#sha256=77586f5deca99bf854dce2bce9e533a90dd97694d190b15bd17c170ef493e2b1 \
  && python3 -m pip install -r requirements.txt
