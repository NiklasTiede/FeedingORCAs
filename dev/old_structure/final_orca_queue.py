#


from subprocess import Popen, PIPE
import subprocess
from pprint import pprint
import os
import re

import time
start_time1 = time.perf_counter()


path = 'testing_orca/many_mols'

# subprocess.call(['rm $(ls -I "*.inp")'], shell=True, cwd=path)
# subprocess.call(['ls', '-la'])


input_output_list = ['ethanol', 'propane', 'water', '2-chloro-2-fluoroethane', '2-methylpropane',
                     'acetaldehyde', 'butane', 'cyclobutane', 'isopropyl_alcohol', 'methanol']
cmds_list = [['/home/niklas/orca/orca', f'{input_output_name}.inp', '>', f'{input_output_name}.out'] for input_output_name in input_output_list]
procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False, cwd=path, text=True) for cmd in cmds_list]

for proc in procs_list:
    print(f'args of this process (subprocess PID={proc.pid}, os PID={os.getpid()}): {proc.args}')
    stdout, stderr = proc.communicate()

    file_name = proc.args[-1][:-4]
    with open(f'{file_name}.xyz', 'r') as f:
        f.readline()

    with open(f'{file_name}.xyz', 'r') as read_file:
        xyz_list = read_file.read().splitlines()

    #pprint(stdout.split('\n')[:50])
    proc.wait()





# --------------------------------------------

# .xyz-file is created and immidiately after creation read out and then deleted (delete other files as well):
import os
os.chdir('testing_orca/many_mols')
with open('2-chloro-2-fluoroethane.xyz', 'r') as read_file:
    lines_list = read_file.read().splitlines()
    for num, line in enumerate(lines_list):
        line1 = line.split()
        if num > 1:
            element, x, y, z = line.split()
            print(num-1, element, float(x), float(y), float(z))

# testset: small molecules and correct calculated coords!





print(f"Script execution finished after {(round(time.perf_counter()-start_time1, 3))} s")



