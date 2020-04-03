set terminal pdf
set output "./data/stat.pdf"
set multiplot layout 3, 1
l_t = 1
l_w = 3

set xlabel "t"
set key top left
p "./data/stat.dat" u 1:2 w l lt l_t lw l_w lc "red" t "vx", "./data/stat.dat" u 1:3 w l lt l_t lw l_w lc "blue" t "vy", "./data/stat.dat" u 1:4 w l lt l_t lw l_w lc "green" t "vz"

set xlabel "t"
set key top left
p "./data/stat.dat" u 1:5 w l lt l_t lw l_w lc "red" t "Ek", "./data/stat.dat" u 1:6 w l lt l_t lw l_w lc "blue" t "Ep", "./data/stat.dat" u 1:7 w l lt 3 lw l_w lc "green" t "Etot"

set xlabel "t"
set key top left
p "./data/stat.dat" u 1:8 w l lt l_t lw l_w lc "red" t "red. temp."

unset multiplot
