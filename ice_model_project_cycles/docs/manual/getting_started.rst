Getting Started
===============

Obtaining ice lens model
------------------------

Get tehpc by cloning it from `gitlab <https://gitlab.ethz.ch>`_::

  > git clone git@gitlab.ethz.ch:youwang/finite-element-simulation-for-coupled-heat-and-moisture-transport.git


Requirements for ice lens model
-------------------------------

The following software are required for tehpc:

- `CMake <https://cmake.org/>`_ (3.1.0 or higher)
  
Optional software for additional features in tehpc:

- `Python3 <https://www.python.org/>`_


Compiling tehpc
-------------------

**Ubuntu-based and macOS systems**:

You can configure and build the project by following these steps after you have the repository on your local machine::

  > mkdir build
  > cd build
  > ccmake ..
  > make

If you would like to run the example as verification, switch on the 'example' option and then make.

  
Running an example
------------------

One example simulation is provided in the `finite-element-simulation-for-coupled-heat-and-moisture-transport/example` folder. To run a simulation, you typically proceed as follows::

  > cd build/example
  > ./example vertical_wall results
  
You may visualize the results with the provided script::

  > python3 results.py

The underlying method and computational options are described in details in :doc:`./user_guide`.
