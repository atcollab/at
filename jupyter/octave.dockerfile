FROM jupyter/base-notebook AS base

USER root
RUN apt-get update -q \
	&& DEBIAN_FRONTEND=noninteractive \
	apt-get install -qy --no-install-recommends \
	git wget gnuplot-x11 openjdk-11-jdk-headless \
	make gcc g++ gfortran \
	gperf flex bison texinfo ghostscript icoutils \
	file libopenblas-dev liblapack-dev librsvg2-bin \
	libqhull-dev libpcre3-dev libncurses-dev libreadline-dev \
	zlib1g-dev bzip2 libhdf5-openmpi-dev libfftw3-dev libglpk-dev \
	libcurl4-openssl-dev libgraphicsmagick++1-dev libfltk1.3-dev \
	libqrupdate-dev libamd2 libsundials-dev libarpack2-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN cd /tmp \
	&& wget https://ftpmirror.gnu.org/octave/octave-6.1.0.tar.gz -O octave.tar.gz \
	&& tar -xvf octave.tar.gz \
	&& cd /tmp/octave-6.1.0 \
	&& ./configure --prefix=/opt/octave CPPFLAGS=-I/usr/include/suitesparse \
	&& cd /tmp/octave-6.1.0 \
	&& make -j$(nproc) \
	&& make install \
	&& cd / \
	&& rm -rf /tmp/octave-6.1.0

ENV PATH="/opt/octave/bin:${PATH}"

FROM base AS at

# TODO: after merge this should clone master at https://github.com/atcollab/at.git
RUN git clone --depth 1 -b octave_compatibilty https://github.com/atcollab/at.git /opt/at \
	&& cd /opt/at/atoctave \
	&& cp -r /opt/at/jupyter/demos $HOME/ \
	&& octave --eval 'bootstrap'

USER jovyan
RUN pip install octave_kernel \
	&& export OCTAVE_EXECUTABLE=$(which octave)

# configure path
RUN cd /opt/at/atoctave && octave --eval 'bootstrap;savepath'
