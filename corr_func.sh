#!/bin/bash

make;
./exec/corr_func;
printf 'set terminal pdf\n set output "./graph/corr_func.pdf"\n p "./data/corr_func.dat" u 1:2 w l t "corr_func"\n exit' | gnuplot
evince ./graph/corr_func.pdf &
