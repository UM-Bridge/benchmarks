================
Tutorial
================

Software setup
========================

The following tutorials assume working installations of:

* Python
* Docker (available from `docker.com <https://www.docker.com/>`_)

However, most steps can be performed using any other programming language with UM-Bridge support.

1: First steps
========================

Minimal model
------------------------

First, install the Python module for UM-Bridge support through PyPI::

    pip install umbridge

You find a minimal UM-Bridge model server written in Python at `UM-Bridge repository <https://github.com/UM-Bridge/umbridge/tree/main/models/testmodel-python/>`_. Download this example server (by git cloning the repository, or just downloading the file itself). Launch it on your machine via::

    python minimal-server.py

The model server is now up and running, waiting to be called by any UM-Bridge client. In particular, it provides a model called "forward" that takes a single 1D vector as input, multiplies it by two, and returns it.


Requesting a model evaluation
------------------------

For now, the server is just idle and waiting. Next, we are going to connect to the model server in order to request a model evaluation.

A basic UM-Bridge client written in Python is ``umbridge-client.py`` from the `UM-Bridge repository <https://www.github.com/UM-Bridge/umbridge/tree/main/clients/python/>`_. Launch it on your machine via::

    python umbridge-client.py http://localhost:4242

``localhost`` in the URL points client to the model running on the same machine. This will lead to an error, indicating that the client is passing an input of wrong dimension to the server.

Take a closer look at the contents of ``umbridge-client.py``. You find an explanation of the Python interface in the `clients section <https://um-bridge-benchmarks.readthedocs.io/en/docs/umbridge/clients.html>`_. Change the input parameter to the dimension the model expects, for example ``param = [[17.4]]``. Running the client again will now yield the expected output from the server, namely ``[[34.8]]``.

Try and request model evaluations for different input parameters; the output should change accordingly.

Note that the client in fact requests the model evaluation twice, the second time passing in a configuration. This model-specific configuration is ignored by the minimal server example, but will become relevant later.


Changing the model
------------------------

Now take a closer look at ``minimal-server.py``. Refer to the `models section <https://um-bridge-benchmarks.readthedocs.io/en/docs/umbridge/models.html>`_ for an explanation of how UM-Bridge models are defined in Python. Play around with the minimal model. For example, you could replace the multiplication by a more interesting operation, or change the model to have a different input and output dimension. Each time, restart the model server and call the modified model from your client to make sure changes take effect as you intend.


Switching out clients
------------------------

Since UM-Bridge is language agnostic, any UM-Bridge client can connect to your minimal model. The syntax is largely the same for any supported language. If you like, follow the instructions in the `clients section <https://um-bridge-benchmarks.readthedocs.io/en/docs/umbridge/clients.html>`_ to call your model from a client in a different language, e.g. C++ or R.

Further, we can use one of the UM-Bridge integrations for specific UQ frameworks. These allow UM-Bridge models to be used seamlessly in each framework by fully embedding UM-Bridge clients (or servers) in the respective architecture. For example, install QMCPy from PyPI::

    pip install qmcpy

You can then run the QMCPy example client from the `UM-Bridge repository <https://www.github.com/UM-Bridge/umbridge/tree/main/clients/python/>`_::

    python qmcpy-client.py http://localhost:4242

It will connect to the running model server as before, and perform a forward UQ solve via Quasi-Monte Carlo. Simply put, it will apply the model to (cleverly chosen) samples from a distribution specified in the client, and output information about the resulting distribution of model outputs. Due to tight integration, this code looks like any other basic QMCPy example; however, it can immediately connect to any (arbitrarily complex) UM-Bridge model. We will take a closer look at UQ frameworks supporting UM-Bridge later.


2: Model containers
========================

Running a pre-defined model container
------------------------

The UM-Bridge benchmark library provides a number of ready-to-run models and UQ benchmark problems. Have a look at the tsunami model; it is part of this documentation site.

Each model in the library is provided as a publicly available Docker image. The Docker image ships not only the model server itself, but also all dependencies and data files it needs. No model-specific setup is required.

As the tsunami model's documentation indicates, it is enough to run the following command to download and run its Docker image::

    docker run -it -p 4242:4242 linusseelinger/model-exahype-tsunami

The model server is now up and running inside a container, waiting to be called by any UM-Bridge client.

Note that only one model server may be running at a given port. So, if you see an error indicating port 4242 is already in use, shut down the existing model server first.

Refer to the tsunami model's documentation again to see what models the model server provides (there may be multiple), and what their properties are. In this case it is a model called ``forward``. This particular model takes a single 2D vector as input, defined to be the location of the tsunami source. It then solves a hyperbolic partial differential equation (PDE) to compute the tsunami propagation. Finally, it returns a single 4D vector containing the main tsunami wave's arrival time and maximum water height at two different locations. This model does not provide any derivatives.


Requesting a model evaluation
------------------------

As before, you can use the minimal Python client to connect to the model (or any other UM-Bridge client).

Apart from input parameters, the client may also choose different configuration options. These are model specific and listed in the respective model's documentation page. For example, the tsunami model allows you to select a finer discretization level by passing ``{"level": 1}`` as configuration. Again, the client documentation gives an example. Be aware that level 2 may take quite long to run.


Switching out models
------------------------

You can use the exact same client as before on any other UM-Bridge model, regardless of model specifics like choice of programming language, build systems, etc. Stop the tsunami model (e.g. via Ctrl + C in its terminal). Instead, run a simple beam benchmark problem::

    docker run -it -p 4243:4243 linusseelinger/benchmark-muq-beam-propagation:latest

This Euler-Bernoulli beam has different input and output dimensions than the tsunami. Running the client again will yield an according error::

    python umbridge-client.py http://localhost:4243

Change the input parameter dimension in ``umbridge-client.py`` to match the new model (e.g. ``param = [[1.02,1.04,1.03]]``), and you should receive its output.

In contrast to the more costly tsunami model, the beam model is fast enough to quickly solve a forward UQ problem on it via QMCPy::

    python3 qmcpy-client.py http://localhost:4243


Accessing model output files
------------------------

Some models may output files in addition to the response the client receives; this is particularly helpful for model debugging. According to its documentation, the tsunami model will write VTK output to the ``/output`` directory if we pass ``{"vtk_output": True}`` as config option.

When launching the model, you can map this directory inside the container to ``~/tsunami_output`` on your machine::

    docker run -it -p 4242:4242 -v ~/tsunami_output:/output linusseelinger/model-exahype-tsunami

After requesting a model evaluation from the client and passing the config option, you can view the output files per time step in your home directory under ``~/tsunami_output`` using paraview or any other VTK compatible visualization tool.


3: Solving UQ problems
========================


Uncertainty propagation
------------------------

We have already looked at uncertainty propagation in passing. Such benchmark problems are essentially equivalent to forward models; however, their documentation specifies a distribution of input parameters, and the goal is to determine (properties of) the resulting distribution of model outputs.

For example, the already mentioned Euler-Bernoulli beam propagation benchmark defines a uniform distribution in three dimesions to sample from. Start the model server now::

    docker run -it -p 4243:4243 linusseelinger/benchmark-muq-beam-propagation:latest

The QMCPy client is already set up to solve the UQ problem defined in the beam benchmark's documentation. Simply run it via::

    python3 qmcpy-client.py http://localhost:4243

Have a closer look at ``qmcpy-client.py``. Try and change the distribution to a different one, e.g. a normal distribution with similar variance. Refer to `QMCPy's documentation <https://qmcpy.readthedocs.io/en/latest/>`_ for details.

Bayesian inverse problems
------------------------

All Bayesian inference benchmarks in the library provide a model named ``posterior`` that maps a model parameter to the log of a Bayesian posterior.
In contrast to propagation benchmarks, the task is to find (properties of) the posterior distribution while only accessing the posterior, and thereby the model, a finite amount of times.
Spin up such a benchmark problem::

    docker run -it -p 4243:4243 linusseelinger/benchmark-analytic-gaussian-mixture

PyMC is a popular package with support for Bayesian inference. It is available via PyPI::

    pip install pymc

The UM-Bridge repository contains a PyMC example client, which you can run as follows::

    python3 pymc-client.py http://localhost:4243

The example uses PyMC's Markov chain Monte Carlo (MCMC) support in order to generate samples from the posterior distribution, only making a finite number of calls to the posterior model. MCMC will explore the parameter space, tending to reject low-posterior samples and accept high-posterior ones. The resulting chain has the posterior distribution as its stationary distribution. Samples from the chain are therefore (correlated) samples from the desired posterior distribution and they may be used to estimate properies of the posterior; the more samples you take, the better the approximation.

This client could also connect to your own model, assuming it provides a model ``posterior`` and has a single 1D output vector (namely the log of the posterior).
The example makes use of PyMC's NUTS sampler to draw samples from the posterior distribution, which is a particular MCMC variant. While this sampler is very efficient, it assumes access to the posterior's gradient. Your model therefore has to provide a gradient implementation for the example to run. Alternatively, you could
switch PyMC to use a different sampler. Refer to `PyMC's documentation <https://www.pymc.io/>`_ for details.


4: Build custom model containers
========================

The easiest way to build your own UM-Bridge model is to create a custom docker container for you model. Docker allows you to package applications, their dependencies, configuration files and/or data to run on Linux, Windows or MacOS systems. They can only communicate with each other through certain channels, we will see more on this later.
In order to create such a docker container you write a set of instructions for building your application. This set of instructions is called a Dockerfile.

Dockerfile structure
------------------------
Writing a Dockerfile is very similar to writing a bash script to build your application. The main advantage is that the Dockerfile will be operating system independant. The main difference is that docker uses certain keywords at the start of each line to denote what type of command you are using.

Before writing our own Dockerfile let's have a look at the Dockerfiles for the two applications we have used in previous steps of the tutorial. The beam propagation benchmark does not have a lot of dependencies. It's Dockerfile can be found `here <https://github.com/UM-Bridge/benchmarks/tree/main/models/muq-beam>`_ .

In addition to the Dockerfile itself the folder contains python files for the applicaton (BeamModel.py and GenerateObservations.py), additional data (ProblemDefinition.h5) and a README. 
We are mainly interested in the Dockerfile itself so let's open it and walk through the components one by one.

On the first line we have::
    
    FROM mparno/muq:latest
    
Here FROM is a keyword we use to define a base image for our application. In this case the model is built on top of the MUQ docker image. The last part `:latest` specifies which version of the container to use.

Next we have::

    COPY . /server

Here COPY is a keyword that specifies we need to copy the server in to the Docker container.

Then we set::

    USER root

The USER keyword can be used to specify which user should be running commands. By default this is root.

Now we need to install any dependencies our application has. In this applications all dependencies can be install using apt and we run::

    RUN apt update && apt install -y python3-aiohttp python3-requests python3-numpy python3-h5py
    
The RUN keyword specifies that the corresponding lines should be executed.

Now we switch user with `USER muq-user` and set the working directory with::

    WORKDIR /server
    
The WORKDIR keyword sets the directory from which all subsequent commands are run. Paths will begin in this directory. If the WORKDIR is not set then `/` is used.

Finally, we run the actual model with::

    CMD python3 BeamModel.py
    
The CMD keyword is also used to execute commands, however, it differs from RUN in that the command is run once container is live. The setup and installation of your application should take place when building the container (use RUN) and the actual model runs should take place once the container is running (use CMD or call this from the umbridge server).

You can also have a look at the Dockerfile for the ExaHyPE tsunami model, which you can find here: `here <https://github.com/UM-Bridge/benchmarks/tree/main/models/exahype-tsunami>`_. This application has more dependencies, and as such a considerably longer Dockerfile, but follows the same steps to install those dependencies one by one. In addition to the keywords described above, this Dockerfile sets environment variables by the `ENV` keyword.

You may notice that this model builds on a base image called `mpioperator/openmpi-builder`. This base image allows you to run MPI commands across docker containers. You can find additional information on this base image `here <https://github.com/kubeflow/mpi-operator>`_.

Comments can be added to a Dockerfile by prepending a `#` character.

Writing your own Dockerfile
------------------------

In order to write your own Dockerfile let's start from the following minimal example.::

    FROM ubuntu:latest

    COPY . /server

    RUN apt update && apt install -y python3-pip

    RUN pip3 install umbridge numpy scipy

    CMD python3 /server/server.py
    
This minimal example assumes a model server is available. Use the model server that you have built in the first part of the tutorial.

Add a file called Dockerfile to your directory. Note that the filename has no extension and is capitalised.

Your Dockerfile should start by building on a base image. As a very basic starting point use ubuntu as your base image::

    FROM ubuntu:latest
    
Alternatively use any other existing image you want to build on.

Next copy the server. Install any standard dependencies your application has::

    RUN apt update && apt install -y python3-pip [your-dependencies]
    
Note:

* python3-pip is needed to install umbridge

* Always remember to run apt update.

* Specify the `-y` option to apt to ensure that apt does not wait for user input.

If you have additional dependencies, add these either by cloning a git repository and installing, or by using the COPY keyword to copy files into your container. 

Install your application. Install umbridge with::
    
    RUN pip3 install umbridge numpy scipy
    
Run the server with::

    CMD python3 /server/server.py.


Building and Running
------------------------

Once you have your Dockerfile you will want to build and run the container. To build the container in your current directory run::

    docker build -t my-model
    
The Dockerfile can also be explicitly set using the -f option. At this stage you may need to go back and modify your Dockerfile because something has gone wrong during the build process.

Once the container is built you can run you model with::

    docker run -it -p 4243:4243 my-model
    
Note that the ports through which your model communicates are specified with the -p option.

It can be useful to check which images currently exist on your computer with::

    docker image ls

Docker images can take up a lot of space and add up quickly. Use `docker image prune` to delete dangling images or `docker image rm` to delete specific images.


(Optional) Uploading to dockerhub
------------------------

Optionally you may want to upload your Dockerfile to dockerhub. This will allow you to build and run by specifying only the name, e.g. ::

    linusseelinger/benchmark-muq-beam-propagation:latest

To push to dockerhub you first need an account. You can set one up at `dockerhub <https://hub.docker.com>`_. Then you can log in on the command line by running::

    docker login
    
Once you are logged in you can push your image to docker hub using::

    docker push my-account/my-model
    
where my-account is your login and `my-account/my-model` is the name of the image you want to push.

