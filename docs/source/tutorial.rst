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

Running a pre-defined model container
------------------------

We begin by running an existing model. The UM-Bridge benchmark library provides a number of ready-to-run models and UQ benchmark problems. Have a look at the tsunami model; it is part of this documentation site.

Each model in the library is provided as a publicly available Docker image. The Docker image ships not only the model server itself, but also all dependencies and data files it needs. No model-specific setup is required.

As the tsunami model's documentation indicates, it is enough to run the following command to download and run its Docker image::

    docker run -it -p 4242:4242 linusseelinger/model-exahype-tsunami

The model server is now up and running inside a container, waiting to be called by any UM-Bridge client.

Refer to the tsunami model's documentation again to see what models the model server provides (there may be multiple), and what their properties are. In this case it is a model called ``forward``. This particular model takes a single 2D vector as input, defined to be the location of the tsunami source. It then solves a hyperbolic partial differential equation (PDE) to compute the tsunami propagation. Finally, it returns a single 4D vector containing the main tsunami wave's arrival time and maximum water height at two different locations. This model does not provide any derivatives.

Requesting a model evaluation
------------------------

For now, the server is just idle and waiting. Next, we are going to connect to the model server in order to request a model evaluation.

Python's UM-Bridge module is available through PyPI::

    pip install umbridge

A basic UM-Bridge client written in Python is ``umbridge-client.py`` from the `UM-Bridge repository <https://www.github.com/UM-Bridge/umbridge/tree/main/clients/python/>`_. Download this example client (by git cloning the repository, or just downloading the file itself). Launch it on your machine via::

    python umbridge-client.py http://localhost:4242

``localhost`` in the URL points client to the model running on the same machine. In the model's terminal, you will now see the log output of the model computation. Finally, the client will receive the model output and print it to terminal.

Take a closer look at the contents of ``umbridge-client.py``. You find an explanation of the Python interface in the clients section of this documentation page. Try and request model evaluations for different input parameters; the output should change accordingly.

Apart from input parameters, the client may also choose different configuration options. These are model specific and listed in the respective model's documentation page. For example, the tsunami model allows you to select a finer discretization level by passing ``{"level": 1}`` as configuration. Again, the client documentation gives an example. Be aware that level 2 may take quite long to run.

Switching out models
------------------------

You can use the exact same client as before on any other UM-Bridge model, regardless of model specifics like choice of programming language, build systems, etc. Stop the tsunami model (e.g. via Ctrl + C in its terminal). Instead, run a simple beam model::

    docker run -it -p 4243:4243 linusseelinger/benchmark-muq-beam-propagation:latest

This Euler-Bernoulli beam has different input and output dimensions than the tsunami. Running the client again will yield an according error::

    python umbridge-client.py http://localhost:4242

Change the input parameter dimension in ``umbridge-client.py`` to match the new model (e.g. ``param = [[0.001]*31]``), and you should receive its output.


Accessing model output files
------------------------

Some models may output files in addition to the response the client receives; this is particularly helpful for model debugging. According to its documentation, the tsunami model above will write VTK output to the ``/output`` directory if we pass ``{"vtk_output": True}`` as config option.

When launching the model, you can map this directory inside the container to ``~/tsunami_output`` on your machine::

    docker run -it -p 4242:4242 -v ~/tsunami_output:/output linusseelinger/model-exahype-tsunami

After requesting a model evaluation from the client and passing the config option, you can view the output files per time step in your home directory under ``~/tsunami_output`` using paraview or any other VTK compatible visualization tool.

Switching out clients
------------------------

We have seen above that the same client can connect to any model. In reverse, any client can connect to a model as well. The syntax is largely the same for any supported language. If you like, follow the instructions in the clients documentation to call a model from a different client integration, e.g. C++ or R.


2: Custom models
========================

You find a minimal UM-Bridge model written in Python at `UM-Bridge repository <https://github.com/UM-Bridge/umbridge/tree/main/models/testmodel-python/>`_. Launch it on your machine via::

    python minimal-server.py

You can now connect to this model as before. It takes a single 1D vector as input, multiplies it by two, and returns it. Refer to the Models section of the documentation for an explanation of how UM-Bridge models are defined in Python. Play around with the minimal model. For example, you could replace the multiplication by a more interesting operation, or change the model to have a different number of inputs and outputs. Call the modified model from your client to make sure the changes take effect as you intend.

3: Solving UQ problems
========================

We now solve an actual UQ problem.

TODO

* simple monte carlo via scipy/numpy random generator?
* QMCPy example on beam propagation ----> best after Aleksei's talk?

4: Build custom container
========================

TODO

* Explanation of Dockerfile structure
* Start from existing Dockerfile
* Adapt to custom model
* How to construct dockerfile by testing container: docker run ... bash

