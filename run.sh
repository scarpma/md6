#!/bin/bash

RUN_DIR=$1

(./exec/main $RUN_DIR &&
./exec/autodiffusion $RUN_DIR  &&
./exec/corr_func $RUN_DIR  &&
python3 graph/stat.py $RUN_DIR &&
python3 graph/autodiffusion.py $RUN_DIR) &&
(open $RUN_DIR/autodiffusion.pdf; open $RUN_DIR/stat.pdf)
#vmd -e /home/scarpma/md6/graph/view.vmd
