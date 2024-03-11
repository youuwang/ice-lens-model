#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as pl
from uqpylab import sessions, display_util

display_util.load_plt_defaults()
uq_colors = display_util.get_uq_color_order()

mysession = sessions.cloud()
uq = mysession.cli
# Reset the session
mysession.reset()
uq.rng(101,'twister')

## read parameters from input file
# def read_parameters(file_path):
#     parameters = {}

#     with open(file_path, 'r') as file:
#         current_section = None

#         for line in file:
#             line = line.strip()

#             if line.startswith('$section'):
#                 _, current_section = line.split(' ', 1)
#                 current_section = current_section.strip()
#                 continue

#             if line.startswith('$'):
#                 # Extract values associated with the parameter
#                 _, *param_values = line.split()
#                 parameters.setdefault(current_section, []).extend(param_values)

#     return parameters

# file_path = 'input/vertical_wall_0.15_203149.99_-14.07_-5.69.inp'  # Replace with the actual path to your file
# parameters = read_parameters(file_path)

# X = [float(x) for x in parameters['material'][2:6]]

items = os.listdir('results')
sub_folder = [item for item in items if os.path.isdir(os.path.join('results', item))]
parameters = []
spacing = []
delta_tot = []
lens_N = []

for folder in sub_folder:
    parameter_list = folder[14:].split('_')
    parameter_temp = [float(x) for x in parameter_list]
    parameter_temp[2] = 10 ** parameter_temp[2]
    parameter_temp[3] = 10 ** parameter_temp[3]
    parameters.append(parameter_temp)
    
    a = open("results/" + folder + "/ub_lb_t.dat", 'r')
    b = open("results/" + folder + "/t_N.dat", 'r')

    a1 = a.readlines()
    b1 = b.readlines()

    del a
    del b

    line_f = b1[-1].split()
    N_lens = int(line_f[1])

    lens_N.append(N_lens)

    x_ub = []
    x_lb = []
    delta = []
    spacing_temp = []
    for i in range(-N_lens, 0):
        lines = a1[i].split()
        x_ub.append(float(lines[0]))
        x_lb.append(float(lines[1]))
    for j in range(len(x_ub)):
        delta.append(x_lb[j] - x_ub[j])
    delta_tot.append(sum(delta))


# spacing = 1000 * np.array(spacing, ndmin = 2).T
delta_tot = 1000 * np.array(delta_tot, ndmin = 2).T
parameters = np.array(parameters, ndmin = 2)

InputOpts = {
    'Marginals': [
        {
        'Name': 'porosity', # porosity
        'Type': 'Uniform',
        'Parameters': [0.15, 0.25]
        },
        {'Name': 'cohesion', # cohesion
        'Type': 'Uniform',
        'Parameters': [2e5, 3e5]
        },        
        {
        'Name': 'permeability at saturation', # permeability at saturated condition
        'Type': 'Uniform',
        'Parameters': [1e-17, 1e-14]
        },
        {'Name': 'pore size', # pore size
        'Type': 'Uniform',
        'Parameters': [1e-6, 1e-5]
        }]
}

myInput = uq.createInput(InputOpts)


PCEOpts = {
    # Select the 'metamodel tool' in UQ[py]Lab
    'Type': 'Metamodel',
    # Choose the Polynomial Chaos Expansion module
    'MetaType': 'PCE',
    # PCE requires an input to be defined
    'Input': myInput['Name']
}
PCEOpts['Degree'] = np.arange(1, 15).tolist()
PCEOpts['ExpDesign'] = {
    'X': parameters.tolist(),
    'Y': np.log(delta_tot).tolist()
}
PCEOpts["Method"] = "OMP"
# Calculate the PCE coefficients
myPCE = uq.createModel(PCEOpts)
# uq.print(myPCE)


# # histogram plot
# fig = pl.figure()
# ax = fig.add_subplot(111)
# y_PCE = uq.evalModel(myPCE, parameters)
# y_PCE = np.exp(y_PCE)
# pl.hist(y_PCE, 20, color = uq_colors[0], alpha = 0.8)
# pl.hist(lens_N, 20, color = uq_colors[1], alpha = 0.8)
# legend_text = ['PCE prediction', 'True model response']
# pl.legend(legend_text, frameon=False, loc="best", fontsize=20)
# ax.set_xlabel('Number of ice lenses', fontsize=20)
# ax.set_ylabel('Counts', fontsize=20)
# ax.tick_params(axis='both', labelsize=20)
# pl.ylim(0, 100)
# # pl.ylim(0, 40)
# pl.show()

# Y-Y plot
y_PCE = uq.evalModel(myPCE, parameters)
y_PCE = np.exp(y_PCE)
pl.scatter(delta_tot, y_PCE, marker = 'o', color = uq_colors[0], alpha = 0.8, s = 10)
pl.plot([0, 250], [0, 250], color = uq_colors[1], alpha = 0.8)
pl.xlabel('$Y_{true}$',fontsize=20)
pl.ylabel('$Y_{PC}$',fontsize=20)
pl.tick_params(axis='both', which='both', labelsize=20)
pl.text(0, 250, 'Total thickness of ice lenses (mm)', fontsize=20, color='black')
pl.grid(False)
pl.show()

# # Sensitivity analysis (Sobol indices)
# PCESobol = {
#     'Type': 'Sensitivity',
#     'Method': 'Sobol',
#     'Sobol': {
#     'Order': 1
#     }
# }

# PCESobolAnalysis = uq.createAnalysis(PCESobol)

# uq.display(PCESobolAnalysis)
# uq.print(PCESobolAnalysis)

