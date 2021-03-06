#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct vec3d {float x, y, z;};
typedef struct vec3d vec;

int t, npartx, nparty, nlayers, npart, write_jump, timesteps, newc, nrun, garb;
float dt, dtsquare, dtdouble, eps, sigma, mu, var, m, a_lattice, last_durata_totale;
float r_max, pot_trunc_perc, BOXL, r_max_squared, shift, penergy, rij2;
float kenergy, simforceij, reduced_temperature, reduced_density;
FILE *output, *cond_in, *cond_in_vel, *logfile, *param, *file_durata_totale;
FILE *inputr, *inputv, *stat, *output_in_cond, *output_in_cond_vel;
vec rij, sumv, vi, rni;

float GaussianNoise(float mu, float var);
void fcc(void);
void load_r(vec *r);
void load_r_2(FILE *inputhere, vec *r);
void write_stat(void);
void write_r(FILE *output_file, vec *r);
void write_vi(FILE *output_file);
void write_durata_totale(void);

float potenergy(float r, float eps, float sigma);
float simforce(float r, float eps, float sigma);
void verlet_periodic(vec *r, vec *ro, vec *a);
void verlet_periodic_write(vec *r, vec *ro, vec *a);
void verlet_periodic_last(vec *r, vec *ro, vec *a);
void compute_forces(vec *r, vec *a);
void compute_forces_stat(vec *r, vec *a);
void eulero(vec *r, vec *ro, vec *a);

