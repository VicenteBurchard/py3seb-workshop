#!/usr/bin/bash -l
source activate py3seb-workshop
jupyter lab --NotebookApp.token=jupyterlab-py3seb-workshop --NotebookApp.allow_origin='http://localhost:3000' & 
myst start
