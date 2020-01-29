'''
This script should only be used if a change 
modified the way netlists are saved.
'''

import os
import shutil

folder = os.path.join(\
    os.path.dirname(__file__),\
    "gui_testing_files")

for subdir, dirs, files in os.walk(folder):
    try:
        shutil.copyfile(
            os.path.join(subdir,'final_after_events_netlist.txt'),
            os.path.join(subdir,'final_netlist.txt')
        )
    except FileNotFoundError:
        pass