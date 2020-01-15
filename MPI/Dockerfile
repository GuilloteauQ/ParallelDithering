FROM ubuntu

MAINTAINER Quentin

RUN apt-get update
RUN apt-get install -qq -y gcc=4:7.4.0-1ubuntu2.3
RUN apt-get install -qq -y openmpi-bin=2.1.1-8
RUN apt-get install -qq -y openmpi-common=2.1.1-8
RUN apt-get install -qq -y openssh-client=1:7.6p1-4ubuntu0.3
RUN apt-get install -qq -y openssh-server=1:7.6p1-4ubuntu0.3
RUN apt-get install -qq -y libopenmpi-dev=2.1.1-8
RUN apt-get install -qq -y make=4.1-9.1ubuntu1
RUN apt-get install -qq -y git

RUN git clone https://github.com/GuilloteauQ/ParallelDithering

WORKDIR /ParallelDithering/MPI

RUN make

ENTRYPOINT ["./run_for_docker.sh"]
