from os.path import dirname, join
import os
import subprocess
import shutil
import sys
repo_directory = dirname(__file__)
docs_directory = join(repo_directory,'docs')
docs_src_directory = join(repo_directory,'docs_src')

# To make sure the make file is found
# by the script
os.chdir(docs_src_directory)

# pip install qucat from the current repository
from pip._internal import main as pipmain
pipmain(['install', repo_directory])

# This is the website url, which is the
# only thing we want to keep from the
# docs directory, we copy it to the docs_src folder
docs_CNAME = join(docs_directory,'CNAME')
docs_src_CNAME = join(docs_src_directory,'CNAME')
try:
    shutil.copy(docs_CNAME,docs_src_CNAME)
except FileNotFoundError:
    print("No CNAME file in docs")

# generate the html code
# the made files are then stored in build/html
os.system('make html')

# replace previous docs with the
# newly built docs
try:
    shutil.rmtree(docs_directory)
except FileNotFoundError:
    print("No docs folder found")
shutil.copytree(join(docs_src_directory,join('build','html')),
    docs_directory)

# Add the website name to the docs folder
shutil.copy(docs_src_CNAME,docs_CNAME)