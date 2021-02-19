FROM debian:bullseye-backports

USER root
RUN apt-get update -q \
	&& DEBIAN_FRONTEND=noninteractive \
	apt install -qy --no-install-recommends \
	git wget build-essential fonts-liberation \
	libpython3-dev python3-pip zip unzip \
	&& DEBIAN_FRONTEND=noninteractive \
	apt install -qy \
	--install-recommends \
	octave liboctave-dev octave-optim octave-parallel \
	&& apt autoremove \
	&& apt autoclean \
	&& rm -rf /var/lib/apt/lists/*

ARG USERNAME
RUN useradd --create-home --shell /bin/bash $USERNAME

USER $USERNAME

ENV PATH="/home/$USERNAME/.local/bin:$PATH" \
	LC_ALL=C.UTF-8 \
	SHELL=/bin/bash

ARG JUPYTERHUB_VERSION
RUN pip3 install jedi==0.17.2 \
	notebook==6.2.0 \
	jupyterlab==3.0.6 \
	jupyterhub==$JUPYTERHUB_VERSION \
	&& jupyter notebook --generate-config \
	&& jupyter lab --generate-config

RUN pip3 install octave_kernel \
	&& export OCTAVE_EXECUTABLE=$(which octave)

ARG DOCKER_NOTEBOOK_DIR AT_REPO AT_BRANCH
RUN git clone --depth 1 -b $AT_BRANCH $AT_REPO /home/$USERNAME/at \
	&& mkdir -p $DOCKER_NOTEBOOK_DIR \
	&& cp -r /home/$USERNAME/at/jupyter/demos $DOCKER_NOTEBOOK_DIR \
	&& cd /home/$USERNAME/at/atoctave \
	&& octave --eval 'bootstrap(:,true);savepath'

RUN cd /home/$USERNAME/at/pyat \
	&& pip3 install -r requirements.txt \
	&& pip3 install matplotlib \
	&& python3 setup.py install --user

WORKDIR $DOCKER_NOTEBOOK_DIR
EXPOSE 8888

COPY scripts/start.sh scripts/start-notebook.sh scripts/start-singleuser.sh /usr/local/bin/
CMD start-notebook.sh --ip 0.0.0.0 --no-browser
