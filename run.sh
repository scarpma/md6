#!/bin/bash

STARTTIME=$(date +%s)
./exec/main &&
./exec/autodiffusion &&
python3 graph/stat.py;
python3 graph/autodiffusion.py;
#gnuplot graph/autodiffusion.gnu;


input="./data/durata_totale.dat";
nrun=$(awk 'NR==2' $input);
stat_final_path=$(printf "./salvati/stat%s.pdf" $nrun)
autodiffusion_final_path=$(printf "./salvati/autodiffusion%s.pdf" $nrun)
stat_dat_final_path=$(printf "./salvati/stat%s.dat" $nrun)
autodiffusion_dat_final_path=$(printf "./salvati/autodiffusion%s.dat" $nrun)

cp "./data/stat.pdf" $stat_final_path;
cp "./data/autodiffusion.pdf" $autodiffusion_final_path;
cp "./data/stat.dat" $stat_dat_final_path;
cp "./data/autodiffusion.dat" $autodiffusion_dat_final_path;

ENDTIME=$(date +%s)
echo "run time = $(($ENDTIME-$STARTTIME)) s" >> ./salvati/run.log
printf "run time = $(($ENDTIME-$STARTTIME)) s\n";
open $stat_final_path &
open $autodiffusion_final_path &
#vmd -e /home/scarpma/md6/graph/view.vmd
