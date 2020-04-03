set terminal pdf
set output "./data/autodiffusion.pdf"
l_t = 1
l_w = 3
set xlabel "t"
set ylabel "varx vary varz"
p "./data/autodiffusion.dat" u 1:2 w l lt l_t lw l_w lc "red" t "varx", "./data/autodiffusion.dat" u 1:3 w l lt l_t lw l_w lc "blue" t "vary", "./data/autodiffusion.dat" u 1:4 w l lt l_t lw l_w lc "green" t "varz"
