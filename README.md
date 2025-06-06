## Molecular Dynamics program implemented in plain C
[![CodeFactor](https://www.codefactor.io/repository/github/scarpma/md6/badge)](https://www.codefactor.io/repository/github/scarpma/md6)

The program was developed for master degree (in Physics) exam "Fisica dei Fluidi Complessi e Turbolenza".

Particles interact as a *Lennard-Jonnes liquid*, with a cutoff range.

*Periodic boundary conditions* are implemented in order to simulate a larger number of particles. Initial positions of particles are the vertex of an *FCC lattice*. This is done in order to avoid particles too near with respect with the initial potential energy and integration time step. Moreover, if the init. pot. energy is right, the FCC vertexes are stable points for the liquid particles and a solid is formed.

The dynamics is integrated through the "*Verlet algorithm*" and some data analysis is done in order to monitor physical quantities such as total momentum and energy. Moreover, a script computes the *auto-diffusion coefficient*.

Balancing the density of particles with the periodic box dimensions is a little difficult and should be faced more seriously.

This repository contains the programs:
- `main.c`
- `autodiffusion.c`
- `corr_func.c`

The autodiffusion coefficient $D$ is defined as: 
$$ \left\langle \left( r(t)-r(0) \right)^2 \right\rangle $$

How to use:
First, compile the programs and prepare a directory where the code will run and save its outputs

```bash
> make
> mkdir run_tmp
```

Then, create a file named `run_tmp/param.in` with the following content (do not modify the order!!)

```
npartx=5
nparty=5
nlayers=5
npart=500
write_jump=10
timesteps=10000
dt=0.0005
eps=10.0
sigma=1.0
mu=0.0
var=1.5
m=1.0
a_lattice=1.5
pot_trunc_perc=0.0005
new_in_cond=1
```

> [!NOTE]
>
> `npart` should be `4*nlayers*npartx*nparty`

Finally, execute the code

```bash
> ./exec/main run_tmp # this runs the MD simulation
> ./exec/autodiffusion run_tmp # this computes the autodiffusion
> ./exec/corr_func run_tmp # this computes the correlation function
```

To analyze the results, use the python scripts:
```bash
> python3 graph/autodiffusion.py run_tmp
> python3 graph/stat.py run_tmp
```



