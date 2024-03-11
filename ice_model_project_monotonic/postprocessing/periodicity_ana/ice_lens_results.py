import os
from lib import r_w as rw
from lib import JobRunner as jr
import numpy as np
import math
import time
import shutil

result_dir = "results_periodic"
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

pars = ['${T_out_ini}', '${T_min}']

def ice_lens(par):

    jobs = []

    for s in range(len(par)):
        par1 = []
        for j in range(len(par[0])):
            par1.append(par[s][j])
        input_file = rw.read_file('vertical_wall_template.inp')

        for i in range(len(par1)):
            input_file = input_file.replace(pars[i], '{}'.format(par1[i]))
        
        input_path = input_dir + "/vertical_wall_{:.2f}_{:.2f}.inp".format(par1[0], par1[1])
        rw.write_file(input_path, input_file)
        result_folder = result_dir + input_path[5:-4]
        if not os.path.exists(result_folder):
            os.mkdir(result_folder)
        args_ex = [input_path[:-4], result_dir + input_path[5:-4]]
        jobs.append([ICE_LENS_EXECUTABLE] + args_ex)

    jobrunner = jr.JobRunner(30)
    jobrunner.run(jobs)
    
    return

X = []
n1 = 20
n2 = 20
x1 = np.linspace(0.5, 1, num = n2)
x2 = np.linspace(-5, -1, num = n1)
for i in range(n2):
    for j in range(n1):
        X.append([x1[i], x2[j]])

start = time.time()
ice_lens(X)
end = time.time()
elapse = end - start
print(elapse)

