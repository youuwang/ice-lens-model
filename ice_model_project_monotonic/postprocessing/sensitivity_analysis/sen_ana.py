#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as pl
from brokenaxes import brokenaxes

# delta_tot = [0.272142, 0.294890, 0.273505, 0.737114, 0.324086]
# lens_N = [0.080127, 0.951589, 0.074983, 0.098141, 0.092248]
# b = [r'$\phi$', r'$f_c$', r'$k_{l,sat}$', r'$r_p$', r'R']
# x = 0.3 * np.arange(len(b))
# width = 0.1
# x1 = [y + width for y in x]
# pl.bar(x, lens_N, color = 'm', width = width)

# pl.bar(x1, delta_tot, color = 'darkorange', width = width)
# pl.xlabel('Material parameters', fontsize = 70)
# pl.ylabel('Total Sobol\' indices', fontsize = 70)
# pl.xticks([val + width / 2 for val in x], b)
# pl.tick_params(axis='both', which='both', labelsize=70, pad = 10)
# pl.legend([r'$N$', r'$\Delta_{tot}$'], frameon=False, loc="upper left", fontsize=70)
# pl.ylim(0, 1)
# ax = pl.gca()
# ax.spines['right'].set_linewidth(2)
# ax.spines['top'].set_linewidth(2)
# ax.spines['left'].set_linewidth(2)
# ax.spines['bottom'].set_linewidth(2)

# pl.subplots_adjust(bottom = 0.17)

# pl.show()


from uqpylab import sessions, display_util

display_util.load_plt_defaults()
uq_colors = display_util.get_uq_color_order()

mysession = sessions.cloud()
uq = mysession.cli
# Reset the session
mysession.reset()
uq.rng(101,'twister')

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
    parameter_temp[4] = 10 ** parameter_temp[4]
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
    for i in range(-N_lens, 0):
        lines = a1[i].split()
        x_ub.append(float(lines[0]))
        x_lb.append(float(lines[1]))
    for j in range(len(x_ub)):
        delta.append(x_lb[j] - x_ub[j])
    delta_tot.append(sum(delta))
# print(lens_N)

# spacing = 1000 * np.array(spacing, ndmin = 2).T
delta_tot = 1e4 * np.array(delta_tot, ndmin = 2).T
parameters = np.array(parameters, ndmin = 2)

InputOpts = {
    'Marginals': [
        {
        'Name': 'porosity', # porosity
        'Type': 'Uniform',
        'Parameters': [0.15, 0.25]
        },
        {'Name': 'cohesion', # cohesion
        'Type': 'Lognormal',
        'Moments': [3e6, 5.5e5]
        },        
        {
        'Name': 'permeability at saturation', # permeability at saturated condition
        'Type': 'Uniform',
        'Parameters': [1e-14, 1e-13]
        },
        {'Name': 'pore size', # pore size
        'Type': 'Uniform',
        'Parameters': [1e-6, 1e-5]
        },
        {'Name': 'grain size', # grain size
        'Type': 'Uniform',
        'Parameters': [1e-4, 8e-4]
        }]
}

myInput = uq.createInput(InputOpts)


PCEOpts_N = {
    # Select the 'metamodel tool' in UQ[py]Lab
    'Type': 'Metamodel',
    # Choose the Polynomial Chaos Expansion module
    'MetaType': 'PCE',
    # PCE requires an input to be defined
    'Input': myInput['Name']
}
PCEOpts_N['Degree'] = np.arange(1, 15).tolist()
PCEOpts_N['ExpDesign'] = {
    'X': parameters.tolist(),
    'Y': np.log(lens_N).tolist()
}
PCEOpts_N["Method"] = "OMP"
# Calculate the PCE coefficients
myPCE_N = uq.createModel(PCEOpts_N)


PCEOpts_delta = {
    # Select the 'metamodel tool' in UQ[py]Lab
    'Type': 'Metamodel',
    # Choose the Polynomial Chaos Expansion module
    'MetaType': 'PCE',
    # PCE requires an input to be defined
    'Input': myInput['Name']
}
PCEOpts_delta['Degree'] = np.arange(1, 15).tolist()
PCEOpts_delta['ExpDesign'] = {
    'X': parameters.tolist(),
    'Y': np.log(delta_tot).tolist()
}
PCEOpts_delta["Method"] = "OMP"
# Calculate the PCE coefficients
myPCE_delta = uq.createModel(PCEOpts_delta)

# # histogram plot
# y_PCE_N = uq.evalModel(myPCE_N, parameters)
# y_PCE_delta = uq.evalModel(myPCE_delta, parameters)

# y_PCE_N = np.exp(y_PCE_N)
# y_PCE_delta = np.exp(y_PCE_delta)

# # fig = pl.figure()
# # ax = fig.add_subplot(111)
# # ax.hist(y_PCE_N, 10, color = uq_colors[0], histtype = 'step',linewidth=4, label = r'PCE prediction of $N$')
# # ax.hist(lens_N, 10, color = uq_colors[1], label = r'True model response of $N$')
# # # ax.hist(y_PCE_delta, 20, color = uq_colors[5], histtype = 'step', linewidth=2, label = r'PCE prediction of $\Delta_{tot}$ (mm)')
# # # ax.hist(delta_tot, 20, color = uq_colors[4], label = r'True model response of $\Delta_{tot}$ (mm)')
# # ax.legend(fontsize = 30, frameon = False, loc = 'best')
# # ax.set_xlabel(r'$N$ or $\Delta_{tot}$', fontsize=35, labelpad = 20)
# # ax.set_ylabel('Counts', fontsize=35, labelpad = 20)
# # ax.spines['top'].set_linewidth(1)
# # ax.spines['right'].set_linewidth(1)
# # ax.spines['bottom'].set_linewidth(1)
# # ax.spines['left'].set_linewidth(1)
# # pl.subplots_adjust(top = 0.95, bottom = 0.17)
# # ax.grid(False)

# # pl.show()


# bax = brokenaxes(
#     ylims=[(0, 100), (200, 300)], 
#     xlims = [(0, 4), (7, 9)],
#     hspace=0.2, 
#     wspace = 0.2,
#     despine = False)

# bax.hist(y_PCE_N, 10, color = uq_colors[0], histtype = 'step',linewidth=4, label = r'PCE prediction of $N$')

# bax.hist(lens_N, 10, color = uq_colors[1], label = r'True model response of $N$')

# bax.hist(y_PCE_delta, 20, color = uq_colors[5], histtype = 'step', linewidth=4, label = r'PCE prediction of $\Delta_{tot}$ (mm)')
# bax.hist(delta_tot, 20, color = uq_colors[4], label = r'True model response of $\Delta_{tot}$ (mm)')
# bax.legend(fontsize = 30, frameon = False, loc = 'best')

# bax.set_xlabel(r'$N$ or $\Delta_{tot}$', fontsize=35, labelpad = 35)
# bax.set_ylabel('Counts', fontsize=35, labelpad = 80)

# bax.first_col[0].set_yticks([200, 300])
# bax.first_col[1].set_yticks([0, 50, 100])
# bax.first_col[0].tick_params(labelsize = 35)
# bax.first_col[1].tick_params(labelsize = 35)

# bax.last_row[0].set_xticks([0, 4])
# bax.last_row[1].set_xticks([7, 9])
# bax.last_row[0].tick_params(labelsize = 35)
# bax.last_row[1].tick_params(labelsize = 35)

# for ax in bax.axs:
#     ax.spines['top'].set_linewidth(1)
#     ax.spines['right'].set_linewidth(1)
#     ax.spines['bottom'].set_linewidth(1)
#     ax.spines['left'].set_linewidth(1)

# pl.subplots_adjust(top = 0.95, bottom = 0.17)
# # pl.tight_layout()

# for handle in bax.diag_handles:
#     handle.remove()
# bax.draw_diags()

# for ax in bax.axs:
#     ax.grid(False)


# pl.show()

# Y-Y plot
y_PCE_N = uq.evalModel(myPCE_N, parameters)
y_PCE_N = np.exp(y_PCE_N)
y_PCE_delta = uq.evalModel(myPCE_delta, parameters)
y_PCE_delta = np.exp(y_PCE_delta)

pl.plot([0, 18], [0, 18], color = 'k', linewidth = 2)
pl.scatter(lens_N, y_PCE_N, marker = 'o', color = uq_colors[1], s = 80)

pl.scatter(delta_tot, y_PCE_delta, marker = 'o', color = uq_colors[4], s = 80)

pl.xlabel('$Y_{true}$',fontsize=35)
pl.ylabel('$Y_{PC}$',fontsize=35)
legend_text = [r'$Y_{PC} = Y_{true}$', r'$N$', r'$\Delta_{tot}$ (mm)']
pl.legend(legend_text, frameon=False, loc="best", fontsize=35)
pl.tick_params(axis='both', which='both', labelsize=35)
pl.grid(False)
pl.subplots_adjust(top = 0.95, bottom = 0.25)

ax = pl.gca()
ax.set_xticks([0, 3, 6, 9, 12, 15, 18])
ax.set_yticks([0, 6, 12, 18])
ax.spines['right'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
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







