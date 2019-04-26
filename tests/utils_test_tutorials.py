import traceback
import nbformat # to convert the notebook into an array of cells
from os.path import join, dirname
from os import devnull
import sys

# Run plt.ion() to avoid hanging on plt.show() calls
import matplotlib.pyplot as plt
plt.ion()

def parse_cell(code):
    lines = code.split('\n')

    valid_lines = ""
    for line in lines:
        if len(line.replace(' ',''))>0:
            if not line.replace(' ','')[0]=='%':
                if 'GUI(' in line and line.rfind(')') != -1:
                    line = line[:line.rfind(')')]+',_unittesting=True'+line[line.rfind(')'):]
                valid_lines += line+'\n'
    return valid_lines



def run_notebook(notebook_name):
    tutorials_folder = join(dirname(dirname(__file__)),'docs','source','tutorials')
    nb = nbformat.read(join(tutorials_folder,notebook_name), 4)
    sys.stdout = open(devnull, "w")
    sys.stderr = open(devnull, "w")
    for c in nb['cells']:
        if c['cell_type'] == 'code':
            to_run = parse_cell(c['source'])
            try: 
                exec(to_run,globals())
            except Exception as e:
                tb = "\nThere was an error executing the following notebook code:\n"
                tb += "==============================\n"
                n=1
                for line in to_run.splitlines():
                    tb += "% 4d %s \n" % (n, line)
                    n += 1
                tb += "==============================\n"
                tb += traceback.format_exc()
                sys.stdout = sys.__stdout__
                return tb
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    return dict(locals())