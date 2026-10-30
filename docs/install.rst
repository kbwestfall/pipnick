.. _install:

Installation
============

``pipnick`` requires Python 3.12 or later.  When installed following the steps
below, all of its dependencies should be installed with it.  Its dependencies
currently include ``numpy``, ``scipy``, ``astropy``, ``astroquery``,
``ccdproc``, ``photutils``, ``reproject``, ``scikit-image``, ``matplotlib``, and
``speclite``.  pipnick therefore leverages much of the scientific python stack
and, in particular, Astropy packages.

Build a new environment
-----------------------

We *strongly* recommend installing ``pipnick`` in its own virtual environment.
This can be done using ``venv`` (preferred; see the `venv documentation
<https://docs.python.org/3/library/venv.html>`__), e.g.:

.. code-block:: console

    python3.13 -m venv ~/.virtualenvs/pipnick
    source ~/.virtualenvs/pipnick/bin/activate

or using `conda <https://docs.conda.io/projects/conda/en/latest/index.html>`__
(see `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__),
e.g.:

.. code-block:: console

    conda create -n pipnick python=3.13
    conda activate pipnick

Install using pip
-----------------

pipnick is not yet released on PyPI, so install the latest version directly from
GitHub:

.. code-block:: console

    pip install "git+https://github.com/UCObservatories/pipnick.git"

Checking the installation
-------------------------

The scripts are installed when the package is installed.  To perform a simple
check that the installation was successful, you can run

.. code-block:: console

    pipnick_version --dependencies

This will print the installed versions of pipnick and all its dependencies.


.. _install-api-key:

The astrometry.net API key
--------------------------

:doc:`astrometry` will submit some images to the astrometry.net web service.  To
use this service, you must register for a free account at
https://nova.astrometry.net and copy your API key from your profile to a
location pipnick can find it.

You should treat these API keys as secret, meaning that you should avoid having
them anywhere that is publicly visible.  We recommend adding the key as an
environmental variable in your shell start-up file (e.g., your ``~/.zshrc`` file
on a modern Mac).  ``pipnick`` assumes the key is assigned to the
``ASTROMETRY_NET_API_KEY`` environment variable.  E.g., add the following line
to your shell start-up file:

.. code-block:: console

    export ASTROMETRY_NET_API_KEY=<your key>

You can also use the ``api_key`` entry in the ``[astrometry_net]`` section of
your ``astroquery`` configuration file (``~/.astropy/config/astroquery.cfg``):

.. code-block:: ini

    [astrometry_net]
    api_key = <your key>

In pipnick, this API key is never a command-line argument and it is never
logged; only the astrometric calibration script needs the key.  pipnick tries to
minimize the number of astrometry.net submissions by leveraging existing
astrometric solutions for a few frames to bootstrap solutions to other frames
observed in the same night.


Developer Installation
----------------------

If you're interested in contributing to the code, you should fork the repository
and then clone and install it as follows:

.. code-block:: console

    git clone https://github.com/<your_github_handle>/pipnick.git
    cd pipnick
    pip install -e ".[dev]"

New contributions require that you check that all tests pass and that the
documentation compiles without any warnings; see below.

Running Tests
~~~~~~~~~~~~~

To test the code, run:

.. code-block:: console

    pytest --pyargs pipnick

The tests use synthetic data and never contact astrometry.net or the Gaia
archive.

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the docs,

.. code-block:: console

    cd docs
    make clean ; make html

The pages are written to ``docs/_build/html``.
