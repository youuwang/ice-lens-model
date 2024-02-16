User's Guide
=================

Finite Difference Method
------------------------

Finite difference method is widely used to solve parabolic partial differential equation like heat conduction equation as shown below.

.. math::

	(\rho c)_{eff} \frac{\partial T}{\partial t} = \frac{\partial (\lambda \frac{\partial T}{\partial x})}{\partial x} + s(x, t)

where :math:`s(x, t)` is the heat souce term which is involved when there is phase change. Central difference method is taken here to express the derivative with respect to space. It should be noted that the solver involved in the current library is just used for 1d case and an explicit method is adopted. If a more flexible time step size is expected, implicit method could be utilized with a linear equation system to be solved, which should be implemented in the :code:`FDM_solver` class or its derived classes. Also, it is possible to extend the FDM to 2d or 3d cases. 


Example
-------

A vertical wall structure with thickness of 1m is solved in the example. Uniform initial temperature profile is taken and temperature at one end is suddenly decreased to a certain value and fixed. Evolution of temperature profile and ice lens formation is predicted within the whole structure. 

The following is the main file for the `example`, which includes all central parts of tehpc.

.. literalinclude:: ../../example/example.cpp
   :language: cpp

	      
Input file
----------

The basic example can be adapted by modifying the input file :code:`vertical_Wall.inp`:

   :language: text

The input file has the following syntax. Keywords start with :code:`$`, and are followed by information. A typical input file contains the following information:

- :code:`$mat_prop` followed by the basic material parameters, including particle density, specific heat capacity, porosity, cohesion, permeability at saturation and characteristic pore throat size.

- :code:`$nodes` followed by the number of nodes used in the simulation.

- :code:`$length` followed by the length or depth of the whole domain investigated.

- :code:`$temperature` followed by the uniform initial temperature within the domain.

- :code:`$bc` followed by the boundary condition applied, namely the temperatures imposed on two ends of the structure.


Running example
---------------

The example is run by the following command::

  > ./example vertical_wall results

The first argument (ie. :code:`vertical_wall`) is the name of the simulation. It will automatically look for an input file with the name of the simulation and extension :code:`.inp` (ie., :code:`vertical_wall.inp`). The second argument (ie. :code:`results`) is the location to write the simulation output. Here, it writes the output directly into the :code:`results` folder.
