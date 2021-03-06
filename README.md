## Molecular Dynamics program implemented in plain C

The program was developed for master degree (in Physics) exam "Fisica dei Fluidi Complessi e Turbolenza".

Particles interact as a *Lennard-Jonnes liquid*, with a cutoff range.

*Periodic boundary conditions* are implemented in order to simulate a larger number of particles. Initial positions of particles are the vertex of an *FCC lattice*. This is done in order to avoid particles too near with respect with the initial potential energy and integration time step. Moreover, if the init. pot. energy is right, the FCC vertexes are stable points for the liquid particles and a solid is formed.

The dynamics is integrated through the "*Verlet algorithm*" and some data analysis is done in order to monitor physical quantities such as total momentum and energy. Moreover, a script computes the *auto-diffusion coefficient*.

Balancing the density of particles with the periodic box dimensions is a little difficult and should be faced more seriously.


