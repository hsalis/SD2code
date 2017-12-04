FROM sd2e/base:ubuntu16

RUN apt-get -y update && apt-get upgrade -y && apt-get install -y python2.7 && apt-get install -y python-pip python-dev build-essential autoconf python-setuptools wget
RUN apt-get install -y git

RUN wget http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.1.tar.gz
RUN tar xvfz mvapich2-2.1.tar.gz
RUN cd mvapich2-2.1/
RUN ./configure --enable-fast=all,O3
RUN make -j 4
RUN cd ..

RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install biopython
RUN pip install scipy
RUN pip install mpi4py

RUN git clone https://github.com/hsalis/SD2code.git

