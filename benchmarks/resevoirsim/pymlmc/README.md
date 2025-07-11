## Synopsis

pymlmc is a demo code for multilevel Monte Carlo computations
in Python.

For mathematical details, see 

    @article{giles2015,
    author = {Giles, M. B.},
    title = {Multilevel {Monte Carlo} methods},
    journal = {Acta Numerica},
    volume = {24},
    month = {5},
    year = {2015},
    pages = {259--328},
    doi = {10.1017/S096249291500001X},
    }

## Authors

P. E. Farrell <patrick.farrell@maths.ox.ac.uk>

M. B. Giles   <mike.giles@maths.ox.ac.uk>

M. Croci      <matteo.croci@maths.ox.ac.uk>

T. Roy        <thomas.roy@maths.ox.ac.uk>

C. Beentjes   <casper.beentjes@maths.ox.ac.uk>

Please contact M. Croci for issues about the code.

## Installation

To download:

    $ git clone https://bitbucket.org/pefarrell/pymlmc.git
    $ cd pymlmc

To run out of the source tree (the easiest way):

    $ export PYTHONPATH=/path/to/pymlmc:$PYTHONPATH

or to install to your user-level Python installation:

    $ python setup.py install --user

or to install to a particular directory:

    $ python setup.py install --prefix=/path/to/directory

## Example

    $ python examples/opre/opre.py

## License

This code is released under the GNU GPL. Anyone requiring a more permissive license
should contact P. E. Farrell and M. B. Giles.