# md6

## Very simple molecular dynamics program implemented in plain C
[![CodeFactor](https://www.codefactor.io/repository/github/scarpma/md6/badge)](https://www.codefactor.io/repository/github/scarpma/md6)

![immagine](https://github.com/user-attachments/assets/b2d00c49-3800-47b1-bf0d-84868c962a5f) ![immagine](https://github.com/user-attachments/assets/1fb2814e-63a5-4b95-9f38-d88501b991bf)

Particles interact as a *Lennard-Jonnes liquid*, with a cutoff range:
```math
u(r_{ij}) = 4\epsilon\left[ \left(\frac{\sigma}{r_{ij}}\right)^{12}-\left(\frac{\sigma}{r_{ij}}\right)^{6} \right]
```
where $u$ is the Lennard-Jonnes potential and $r_{ij}=\left|x_i-x_j\right|$ is the distance between particles. For $r_{ij}>r_{cut}$, $u = 0$. Particle positions are evolved solving Newtonian mechanics:
```math
F_i=m\frac{d^2r}{dt^2}
```
```math
F_i=-\nabla_i u(r).
```

These equations are numerically solved using the **Verlet algorithm**:
```math
r_i(t+\Delta t)=2r_i(t)-r_i(t-\Delta t)+\frac{F_i(t)}{m}\Delta t^2.
```

**Periodic boundary conditions** are implemented in order to simulate a larger number of particles. Initial positions of particles are the vertex of an **FCC lattice**. This is done in order to avoid particles too near with respect with the initial potential energy and integration time step. Moreover, if the init. pot. energy is right, the FCC vertexes are stable points for the liquid particles and a solid is formed.

> [!NOTE]
> Some data analysis is performed in order to monitor physical quantities such as total momentum and energy (`<run_dir>/stat.dat`, `graph/stat.py` and `<run_dir>/stat.pdf`). Please, consider that given the cutoff range $r_{cut}$ and the periodic conditions, energy and momentum may bot be measured / conserved appropriately. **Consider this code as a toy model**. 

The program computes the **auto-diffusion coefficient** (`autodiffusion.c`, `graph/autodiffusion.py` and `<run_dir>/autodiffusion.pdf`).

> [!NOTE]
> Balancing the density of particles with the periodic box dimensions is a little difficult and should be faced more seriously. Take this as a toy example.

The autodiffusion coefficient $D$ is defined as: 
```math
\sigma^2_r(t) = \frac{1}{N_p}\sum_{i=1}^{N_p}\left\langle \left( r_i(t)-r_i(0) \right)^2 \right\rangle
```
```math
D = \lim_{t\to\infty}\frac{1}{6t}\sigma^2_r(t)
```
where $N_p$ is the number of particles of the simulation. If $D\simeq0$, it means that particles are almost fixed w.r.t. their initial positions, i.e. we are simulating a solid. If $D>0$, then we are simulating a fluid.

This repository contains the programs:
- `main.c`
- `autodiffusion.c`
- `corr_func.c`

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
```
and open the plotted pdf files in `run_tmp/autodiffusion.pdf` and `run_tmp/stat.pdf`

> [!NOTE]
> The program was developed as exam project for "Fisica dei Fluidi Complessi e Turbolenza"
