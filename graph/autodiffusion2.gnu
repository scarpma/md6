set terminal pdf

set output "autodiffusion.pdf"
p "./data/autodiffusion.dat" u 1:2, "./data/autodiffusion.dat" u 1:3, "./data/autodiffusion.dat" u 1:4
