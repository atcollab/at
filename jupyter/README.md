This folder contains a docker files for running a
[JupyterHub](https://jupyterhub.readthedocs.io/en/stable/) server.
It borrows from [jupyterhub-deploy-docker](https://github.com/jupyterhub/jupyterhub-deploy-docker)
repo and provides multi-user Jupyter Notebook environment for small classes, teams or departments.

Requirements for running server:

- [`docker`](https://www.docker.com/)
- [`docker-compose`](https://github.com/docker/compose)
- [`make`](https://www.gnu.org/software/make/)

## Basic usage
### Start JupyterHub

```bash
make
docker-compose up -d
```

Server will listen on port `443`.
You may connect to it form browser using one of the lines below:

- <https://127.0.0.1> -- for local installations
- <https://1.1.1.1> -- where `1.1.1.1` is your real IP. Use this when you have a real IP address but
- <https://my-domain.org> -- where `my-domain.org` is your domain name

Note that browsers will complain if your certificate was not verified by [Certificate
authority](https://en.wikipedia.org/wiki/Certificate_authority).
For a small server like this it's easier to check certificate manually and then click on `Accept the
risk` button in your browser.

### Stop JupyterHub

```bash
docker-compose down
```

You will **not** lose your data is you'll stop the server.

### Inspect logs

```bash
docker-compose logs --follow
```

### Remove saved data

```bash
make purge
```

## First run

At first run you will be prompted to enter basic information for the self signed certificate which
will be used for HTTPS connections.
It will be located in `secrets` folder.

After that docker will pull and build required images, which may take a while.

### Setup admin

The first thing you need to do is to setup admin password. For that proceed to the Signup page and
create a user with Username `admin` (you may change it in admin panel later).
Now you can log into the server under `admin` credentials.

### Create a new user

This setup uses [native authenticatior](https://github.com/jupyterhub/nativeauthenticator) for
JupyterHub.

In order to add a new user:

1. User must proceed to the Signup page and enter credentials;
2. Admin must access <https://127.0.0.1/hub/authorize> page (there's no link to it yet. See
[this issue](https://github.com/jupyterhub/nativeauthenticator/pull/79)) and click `Authorize`
button next to the user name.
3. Now new user may login into the server with provided credentials.
