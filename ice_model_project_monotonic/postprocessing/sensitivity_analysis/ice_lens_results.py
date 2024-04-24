import os
from lib import r_w as rw
from lib import JobRunner as jr
import numpy as np
import math
import time
from uqpylab import sessions
import shutil

mysession = sessions.cloud()
uq = mysession.cli
# Reset the session
mysession.reset()
uq.rng(101,'twister')

result_dir = "results"
input_dir = "input"
# get absolute file path to DELPHIN solver (assumed to be inside this directory)
ICE_LENS_EXECUTABLE = './example'
wd = '/home/youwang/Desktop/C++ project/ice_model_project/sensitivity_analysis'

if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

if os.path.exists(input_dir):
    shutil.rmtree(input_dir)

os.mkdir(result_dir)
os.mkdir(input_dir)

pars = ['${porosity}', '${cohesion}', '${kl_eff}', '${pore_size}', '${grain_size}']

def ice_lens(par):

    jobs = []

    for s in range(len(par)):
        flag = 1
        par1 = []
        for j in range(len(par[0])):
            par1.append(par[s][j])
        input_file = rw.read_file('vertical_wall_template.inp')

        for i in range(len(par1)):
            input_file = input_file.replace(pars[i], '{}'.format(par1[i]))
        
        input_path = input_dir + "/vertical_wall_{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.2f}.inp".format(par1[0], par1[1], math.log10(par1[2]), math.log10(par1[3]), math.log10(par1[4]))
        rw.write_file(input_path, input_file)
        result_folder = result_dir + input_path[5:-4]
        if not os.path.exists(result_folder):
            os.mkdir(result_folder)
        args_ex = [input_path[:-4], result_dir + input_path[5:-4]]
        jobs.append([ICE_LENS_EXECUTABLE] + args_ex)

    jobrunner = jr.JobRunner(30)
    jobrunner.run(jobs)
    
    return

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
        }
        ]
}

myInput = uq.createInput(InputOpts)

X = uq.getSample(N = 400, Method = 'LHS')

start = time.time()
ice_lens(X)
end = time.time()
elapse = end - start
print(elapse)

