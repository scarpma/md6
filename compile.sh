#!/bin/bash

gcc -o ./exec/main main.c verlet_periodic.c in_cond.c -lm;
gcc -o ./exec/autodiffusion autodiffusion.c -lm;
