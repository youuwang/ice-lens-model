#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import griddata

results = "results_periodic"
items = os.listdir(results)
sub_folder = [item for item in items if os.path.isdir(os.path.join(results, item))]
parameters = []
lens_N = []

for folder in sub_folder:
    d = []
    parameter_list = folder[14:].split('_')
    parameter_temp = [float(x) for x in parameter_list]

    parameters.append(parameter_temp)
    
    a = open(results + "/"  + folder + "/t_N.dat", 'r')
    b = open(results + "/"  + folder + "/ub_lb_t.dat", 'r')

    a1 = a.readlines()
    b1 = b.readlines()

    del a
    del b

    line_f = a1[-1].split()
    N_lens = int(line_f[1])

    for i in range(-N_lens, 0):
        lines = b1[i].split()
        x_ub = float(lines[0])
        x_lb = float(lines[1])
        d.append(x_lb - x_ub)

    lens_N.append(N_lens)
parameters = np.array(parameters, ndmin = 2)
lens_N = np.array(lens_N)

grid_x, grid_y = np.mgrid[min(parameters[:, 0]):max(parameters[:, 0]):100j, min(parameters[:, 1]):max(parameters[:, 1]):100j]
grid_colors = griddata(parameters, lens_N, (grid_x, grid_y), method='linear')

mask = np.isnan(grid_colors)
grid_colors[mask] = griddata(parameters, lens_N, (grid_x[mask], grid_y[mask]), method='nearest')
image = pl.imshow(grid_colors.T, extent=(min(parameters[:, 0]), max(parameters[:, 0]), min(parameters[:, 1]), max(parameters[:, 1])), origin='lower', cmap='viridis')
pl.gca().set_aspect('auto', adjustable='box')
cbar = pl.colorbar(label = 'Number of ice lenses')
cbar.ax.tick_params(labelsize=30)
cbar.set_label('Number of ice lenses', fontsize=30)
image.set_clim(1, 70)
ticks = [1, 24, 47, 70]
cbar.set_ticks(ticks)

pl.xlabel(r'$T_0$ ($^\circ \mathrm{C}$)', fontsize = 30)
pl.ylabel(r'$T_1$ ($^\circ \mathrm{C}$)', fontsize = 30)
pl.xlim(0.5, 1)
pl.ylim(-5, -1)

pl.tick_params(axis='both', labelsize = 30)
pl.show()