JOSS
====

.. title:: 1d ice lens model coupled with heat conduction

:Authors: - You Wang, you.wang@ifb.baug.ethz.ch

:Affiliations: - Computational Mechanics for Building Materials (CMBM), Institute for Building Materials, ETHz

Summary
-------

Damages of old structures built with porous materials could be caused by multiple reasons, among which freezing plays a role, and it has been challenging and necessary to get a fundamental understanding of the underlying mechanism for freezing-related fractures in porous medium. Additionally, moisture transport and heat transfer both have great effects on the location of phase changes. Thus, it will be interesting to offer a software to implement such an ice lens model which couples the above physical processes. 

This software intends to predict ice lens formation and growth under 1d case, considering heat conduction process within porous media. Typically, we have a uniform temperature profile initially within the satured material and then one end of the domain is cooled down suddenly. The software could work from cases without any ice to cases with periodic ice lens formation. 


Statement of need
-----------------

Finite Difference Method (FDM) is widely used to solve parabolic partial differential equation like heat conduction equation as shown below using apparent heat capacity [1]_.

.. math::

	C_{a} \frac{\partial T}{\partial t} = \frac{\partial (\lambda \frac{\partial T}{\partial x})}{\partial x}

where :math:`C_{a}` is apparent heat capacity and could be given as

.. math::

	C_{a} = (\rho c)_{eff} - \rho_{i} L \phi \frac{\partial S_{i}}{\partial T}

where :math:`L` is latent heat of fusion for water; :math:`\phi` represents porosity; and :math:`S_{i}` is ice saturation. Explicit and implicit schemes have both been well developed for discretization of the partial differential equation. While explicit method is easier to implement, the implicit one enables a more flexible time step size without constraints for convergency. Right now, the solver in the software uses the explicit method.

For formation of ice lenses in porous media, Rempel et al. (2004) [2]_ proposed a model to predict ice lens formation and growth in a frost heave based on thermodynamics and force equilibrium, considering premelting dynamics between soil particles and ice. Rempel (2007) [3]_ further investigated the Stefan problem where the heat conduction equation is simplified by ignoring differences in material properties.

The above work about ice lens models, however, fails to take into account heat conduction process in a more realistic way with both discontinuities of material properties and heat source due to phase changes considered. For example, in Rempel (2004) [2]_, a linear temperature profile is assumed during the whole process, which is always not the case. Thus, a more realistic ice lens model is required in order to gain a more precise understanding of the fractures formed in porous materials due to freezing, which is the main research purpose of this software.

Acknowledgement
---------------

This project is supported by the course "Towards Efficient and High-Performance Computing in engineering" given by Prof. David Kammer. Sincere thanks to the lecturer and all teaching assistants for the helpful suggestions and feedbacks.

References
----------

.. [1] Tubini, Niccolò; Gruber, Stephan; Rigon, Riccardo, "A method for solving heat transfer with phase change in ice or soil that allows for large time steps while guaranteeing energy conservation," *The Cryosphere*, 15, 2021, pp. 2541–2568.  https://doi.org/10.5194/tc-15-2541-2021.

.. [2] Rempel, Alan; Wettlaufer, J.S.; Worster, M. Grae, "Premelting dynamics in a continuum model of frost heave," *Journal of Fluid Mechanics*, 498, 2004, pp. 227-244. https://doi.org/10.1017/S0022112003006761.

.. [3] Rempel, Alan, "Formation of ice lenses and frost heave," *Journal of Geophysical Research*, 112, 2007, F02S21. https://doi.org/10.1029/2006JF000525. 
