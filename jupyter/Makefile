# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

include .env

.DEFAULT_GOAL=all

.PHONY: network
network:
	@docker network inspect $(DOCKER_NETWORK_NAME) >/dev/null 2>&1 || docker network create $(DOCKER_NETWORK_NAME)

.PHONY: volumes
volumes:
	@docker volume inspect $(DATA_VOLUME_HOST) >/dev/null 2>&1 || docker volume create --name $(DATA_VOLUME_HOST)
	@docker volume inspect $(DB_VOLUME_HOST) >/dev/null 2>&1 || docker volume create --name $(DB_VOLUME_HOST)

secrets/postgres.env:
	@mkdir -p secrets
	@echo "Generating postgres password in $@"
	@echo "POSTGRES_PASSWORD=$(shell openssl rand -hex 32)" > $@

secrets/jupyterhub.key: secrets/jupyterhub.crt

secrets/jupyterhub.crt:
	@mkdir -p secrets
	@openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout \
	secrets/jupyterhub.key -out secrets/jupyterhub.crt

.PHONY: notebook_image
notebook_image: at.dockerfile
	docker build --build-arg USERNAME=$(USERNAME) \
		--build-arg JUPYTERHUB_VERSION=$(JUPYTERHUB_VERSION) \
		--build-arg DOCKER_NOTEBOOK_DIR=$(DOCKER_NOTEBOOK_DIR) \
		--build-arg AT_BRANCH=$(AT_BRANCH) \
		--build-arg AT_REPO=$(AT_REPO) \
		-t $(DOCKER_NOTEBOOK_IMAGE) -f $^ .

.PHONY: build
build: secrets/jupyterhub.crt secrets/jupyterhub.key secrets/postgres.env network volumes
	docker-compose build

.PHONY: all
all: build notebook_image

.PHONY: purge
purge:
	docker-compose down
	-@docker network inspect $(DOCKER_NETWORK_NAME) >/dev/null 2>&1 && docker network rm $(DOCKER_NETWORK_NAME)
	-@docker volume inspect $(DATA_VOLUME_HOST) >/dev/null 2>&1 && docker volume rm $(DATA_VOLUME_HOST)
	-@docker volume inspect $(DB_VOLUME_HOST) >/dev/null 2>&1 && docker volume rm $(DB_VOLUME_HOST)
