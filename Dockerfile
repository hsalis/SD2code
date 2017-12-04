FROM sd2e/base:ubuntu16

RUN apt-get -y update && apt-get upgrade -y && apt-get install -y python2.7 && apt-get install -y python-pip python-dev build-essential autoconf python-setuptools wget
RUN apt-get install -y git

# RUN wget http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.1.tar.gz
# RUN tar xvfz mvapich2-2.1.tar.gz
# RUN cd mvapich2-2.1 && autoconf && ./configure --enable-fast=all,O3 && make -j 4 && cd ..

# Install SD2e CLI tools
RUN curl -L https://raw.githubusercontent.com/sd2e/sd2e-cli/master/install/install.sh | sh
RUN bash ~/.bashrc

RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install biopython
RUN pip install scipy
# RUN pip install mpi4py

# This code will be needed if github repo is private
# Make ssh dir
# RUN mkdir /root/.ssh/

# Copy over private key, and set permissions
# ADD id_rsa /root/.ssh/id_rsa

# Create known_hosts
# RUN touch /root/.ssh/known_hosts
# Add bitbuckets key
# RUN ssh-keyscan github.com >> /root/.ssh/known_hosts

# This github repo is public
# RUN git clone https://github.com/hsalis/SD2code.git

ADD . *
