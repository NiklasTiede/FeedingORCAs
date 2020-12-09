from src.my_timer import time_converter

# TODO: putting input files in a queue for ORCA

# feeding every CPU core with one Molecule and observe the process, time meas. -> create time estimator from those data

# extract certain data from .output files and write them into a database/ my classes instances
# database has features to visualize properties of the molecules ??

# feature: queue for a selectable time (like 22:00) wake up OS from sleep and start: first: adding .py script to win
# scheduler

# parallel computation. works perfectly !
from subprocess import Popen, PIPE
from src.my_timer import time_converter
import time

start_time1 = time.perf_counter()  # total time count

path = '3_molecules'
input_output_list = ['water']
#input_output_list = ['ethanol', 'propane', 'water', '2-chloro-2-fluoroethane', '2-methylpropane',
                     #'acetaldehyde', 'butane', 'cyclobutane', 'isopropyl_alcohol', 'methanol']
cmds_list = [f'orca {input_output_name}.inp > {input_output_name}.out' for input_output_name in input_output_list]
procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, cwd=path) for cmd in cmds_list]
print(procs_list)
for proc in procs_list:
    proc.wait()
    #print(proc.wait())  # 0 -> no errors
    #print(PIPE)  # -> -1

print(f"Script execution finished after {time_converter(round(time.perf_counter()-start_time1, 3))}")

# -------------

# # script for changing the header of a .input file: adressing just line number (here 1 and 2). works perfectly !
# import time
# import os
# print()
#
# start_time1 = time.perf_counter()
#
# os.chdir('writing_orca_commands')
#
# with open("propane.inp", 'r') as read_file:
#     file_as_list = read_file.readlines()
#
# file_comment = '# Basic Mode\n'
# file_command = '! B3LYP def2-SVP Opt\n'
#
# with open("test.inp", 'w') as write_file:
#     file_as_list[0], file_as_list[1] = file_comment, file_command
#     list_back_to_string = ''
#     for line in file_as_list:
#         list_back_to_string += line
#     write_file.write(list_back_to_string)
#
# print(f"\nScript execution finished after {time_converter(round(time.perf_counter()-start_time1, 3))}")
#
# # use: converting calculated .xyz file to new .inp files
# # use: adding commands to other files (maybe smiled converted xyz data)

# ------------
# script for saving data from certain output files. search with key_string
# (data like molecule properties etc?)

import time
import re
import os

#print()
start_time1 = time.perf_counter()

# comment_lines = '# Basic Mode\n'
# keyword_lines = '! B3LYP def2-SVP Opt\n'
# file_block = ''
os.chdir('writing_orca_commands')

with open("propane.inp", 'r') as read_file:
    file_as_list = read_file.readlines()

file_comment = '# Basic Mode\n'
file_command = '! B3LYP def2-SVP Opt\n'

with open("test.inp", 'w') as write_file:
    file_as_list[0], file_as_list[1] = file_comment, file_command
    list_back_to_string = ''
    for line in file_as_list:
        list_back_to_string += line
    write_file.write(list_back_to_string)


def keyword_command_inserter(file_name, path, new_comment, new_keyword):
    ''' Replaces the commands/comment of an .inp file (which contains the '#' keyword
    for comment and '!' for command of the last job) with cartesian coordinates
    by new commands/comment for a new ORCA job.
    example: # file created by avogadro -> # command replacement (file was edited at 18.04.2020)
             ! B3LYP def2-SVP Opt       -> !
    '''
    os.chdir(path)
    with open(f"{file_name}.inp", 'r') as read_file:
        file_as_list = read_file.readlines()
    for line in enumerate(file_as_list):
        if '#' in line[1]:
            file_as_list[line[0]] = f"# {new_comment}\n"  # plz dont use '!'
        if '!' in line[1]:
            file_as_list[line[0]] = f"! {new_keyword}\n"
    with open("test.inp", 'w') as write_file:
        list_back_to_string = ''
        for line in file_as_list:
            list_back_to_string += line
        write_file.write(list_back_to_string)


list1 = ['ab', '!Puk']

for line in enumerate(list1):
    if '!' in line[1]:
        list1[line[0]] = "new command1"  # all lines containing '!' are changed


print(list1)

print(f"\nScript execution finished after {time_converter(round(time.perf_counter()-start_time1, 3))}")

# --------------
# how to convert SMILES files into xyz coords (using pybel!)

# keyword_command_inserter(file_name, 'writing_orca_commands', 'Basic Mode', 'B3LYP def2-SVP Opt')


# def keyword_command_inserter-old_inp(file_name, path, new_comment, new_keyword, new_block=str()):
#     pass
#
#
# def command_inserter-xyz_file(file_name, path, comment, keyword, block):
#     pass


def main():
    """ testing the code. """
    pass


if __name__ == '__main__':
    main()






















