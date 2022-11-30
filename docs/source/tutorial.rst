================
Tutorial
================

Software setup
========================

The following tutorials assume working installations of:

* Python
* Docker (available from `docker.com <https://www.docker.com/>`)

However, most steps can be performed using any other programming language with UM-Bridge support.

1: First steps
========================

Running an existing model container
------------------------

We begin by running an existing model. The UM-Bridge benchmark library provides a number of ready-to-run models and UQ benchmark problems. Have a look at the tsunami model; it is part of this documentation site.

Each model in this library is provided as a publicly available Docker image. The Docker image ships not only the model server itself, but also all dependencies and data files it needs. Therefore, no model-specific setup is required.

As the tsunami model's documentation indicates, it is enough to run the following command to download and run the Docker image::

    docker run -it -p 4242:4242 linusseelinger/model-exahype-tsunami

The model server is now up and running inside a container, waiting to be called by any UM-Bridge client.

Refer to the model's documentation again to see what models the model server provides (there may be multiple), and what their properties are. In this case it is a model called ``forward``. This particular model takes a single 2D vector as input, defined to be the location of the tsunami source. It then solves a hyperbolic partial differential equation (PDE) to compute the tsunami propagation. Finally, it returns a single 4D vector containing the main tsunami wave's arrival time and maximum water height at two different locations. This model does not provide any derivatives.

Requesting a model evaluation
------------------------

For now, the server is just idle and waiting. We are now going to connect to the model server in order to request a model evaluation.

Python's UM-Bridge module is available through PyPI::

    pip install umbridge

A basic UM-Bridge client written in Python is ``umbridge-client.py`` from the `UM-Bridge repository <https://www.github.com/UM-Bridge/umbridge/tree/main/clients/python/>`_. Download this example client (by git cloning the repository, or just downloading the file itself). Launch this example client for Python on your machine via::

    python3 umbridge-client.py http://localhost:4242




    docker run -it -p 4242:4242 -v ~/tsunami_output:/output linusseelinger/model-exahype-tsunami


========================
TODO
========================



* Run existing client from repo on tsunami model (follow quickstart)
* inspect output
* change input parameters
* change config
* Switch out model for other container from benchmark lib
* Optionally try out different client examples (C++, R, curl)

2: Build custom model
========================
* Model class in own program with simple analytic model
* Then move model class into separate program
* Access model server from own client and QMCPy

3: Solving UQ problems
========================
* simple monte carlo via scipy/numpy random generator?
* QMCPy example on beam propagation ----> best after Aleksei's talk?

4: Build custom container
========================
* Explanation of Dockerfile structure
* Start from existing Dockerfile
* Adapt to custom model
* How to construct dockerfile by testing container: docker run ... bash

