FROM bioconductor/bioconductor_docker
USER root

RUN apt-get update -qq && apt-get -y --no-install-recommends install cmake subversion git bc
RUN cd /opt/ && git clone https://github.com/LisaNicvert/camtrapHawkes.git
