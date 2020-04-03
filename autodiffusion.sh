#!/bin/bash

make;
./exec/autodiffusion;
printf 'set terminal pdf\n set output "./graph/autodiffusion.pdf"\n p "./data/autodiffusion.dat" u 1:2 w l t "x^2", "./data/autodiffusion.dat" u 1:3 w l t "y^2", "./data/autodiffusion.dat" u 1:4 w l t "z^2", "./data/autodiffusion.dat" u 1:($2 + $3 + $4) w l t "r^2"\n exit' | gnuplot
evince ./graph/autodiffusion.pdf &
