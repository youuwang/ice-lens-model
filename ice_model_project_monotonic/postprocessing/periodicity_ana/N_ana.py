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
    
    a = open(results + "/"  + folder + "/t_N.dat", 'r')
    b = open(results + "/"  + folder + "/ub_lb_t.dat", 'r')

    a1 = a.readlines()
    b1 = b.readlines()

    del a
    del b

    if len(a1) == 0:
        continue

    parameter_list = folder[14:].split('_')
    parameter_temp = [float(x) for x in parameter_list]

    parameters.append(parameter_temp)

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
cbar = pl.colorbar(orientation='horizontal', label = 'Number of ice lenses', shrink = 0.8, pad = 0.2)
cbar.ax.tick_params(labelsize=70)
cbar.set_label('Number of ice lenses', fontsize=70)
image.set_clim(1, 29)
ticks = [1, 15, 29]
cbar.set_ticks(ticks)

pl.xlabel(r'$T_0$ ($^\circ \mathrm{C}$)', fontsize = 70)
pl.ylabel(r'$T_1$ ($^\circ \mathrm{C}$)', fontsize = 70, labelpad = 5)
pl.xlim(1, 5)
pl.ylim(-10, -1)
pl.gca().set_yticks([-10, -1])
pl.gca().set_xticks([1, 3, 5])

pl.tick_params(axis='both', pad = 21, labelsize = 70)
pl.subplots_adjust(top = 0.95, bottom = 0.1)

pl.show()