# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
FROM debian:testing

WORKDIR /usr

ENV LC_ALL=C.UTF-8

RUN apt-get update -q \
	&& DEBIAN_FRONTEND=noninteractive \
	apt-get install -qy --no-install-recommends \
        build-essential python3-pip python3-psycopg2 \
        libpq-dev npm \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN npm install -g configurable-http-proxy

# Install dockerspawner, authentificators
ARG JUPYTERHUB_VERSION
RUN pip3 install jupyterhub==$JUPYTERHUB_VERSION \
	jupyterlab==3.0.5 \
        jupyterhub-nativeauthenticator==0.0.7 \
        jupyterhub-dummyauthenticator==0.3.1 \
        dockerspawner==0.11.1

# Copy TLS certificate and key
ENV SSL_CERT /srv/jupyterhub/secrets/jupyterhub.crt
ENV SSL_KEY /srv/jupyterhub/secrets/jupyterhub.key
COPY ./secrets/*.crt $SSL_CERT
COPY ./secrets/*.key $SSL_KEY
COPY ./jupyterhub_config.py /srv/jupyterhub/jupyterhub_config.py
RUN chmod 700 /srv/jupyterhub/secrets && \
    chmod 600 /srv/jupyterhub/secrets/*

CMD jupyterhub
