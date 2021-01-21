FROM debian:bullseye-backports

USER root
RUN apt-get update -q \
	&& DEBIAN_FRONTEND=noninteractive \
	apt-get install -qy --no-install-recommends \
	git wget build-essential fonts-liberation \
	python3-pip jupyter-notebook \
	&& DEBIAN_FRONTEND=noninteractive \
	apt-get install -qy --install-recommends \
	octave liboctave-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

ENV LC_ALL=C.UTF-8

RUN pip3 install octave_kernel \
	&& export OCTAVE_EXECUTABLE=$(which octave)

ARG USER="atuser"

RUN useradd --create-home --shell /bin/bash $USER

USER $USER

# TODO: after merge this should clone master at https://github.com/atcollab/at.git
RUN git clone --depth 1 -b octave_compatibilty https://github.com/atcollab/at.git /home/$USER/at \
	&& mkdir /home/$USER/workdir \
	&& cp -r /home/$USER/at/jupyter/demos /home/$USER/workdir \
	&& cd /home/$USER/at/atoctave \
	&& octave --eval 'bootstrap;savepath'

WORKDIR /home/$USER/workdir
EXPOSE 8888
CMD jupyter notebook \
	--no-browser \
	--ip=0.0.0.0 \
	--port=8888 \
	--port-retries=0 \
	--notebook-dir=/home/atuser/workdir
