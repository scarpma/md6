#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// useful for debbugging
#define LOG_LOCATION() printf("File: %s, Line: %d\n", __FILE__, __LINE__)

struct vec3d {float x, y, z;};
typedef struct vec3d vec;

int reproducible;
int t, npartx, nparty, nlayers, npart, write_jump, timesteps, newc, nrun, garb;
float dt, dtsquare, dtdouble, eps, sigma, mu, var, m, a_lattice, last_durata_totale;
float r_max, pot_trunc_perc, BOXL, r_max_squared, shift, penergy, rij2;
float kenergy, simforceij, reduced_temperature, reduced_density;
#define STRLEN 256
char paramfilepath[STRLEN], logfilepath[STRLEN];
char statfilepath[STRLEN], statfilepath[STRLEN];
char coordfilepath[STRLEN], restartcoordfilepath[STRLEN];
char restartvelfilepath[STRLEN], autodiffusionfilepath[STRLEN];
char totdurationfilepath[STRLEN], corrfuncfilepath[STRLEN];

FILE *paramfile, *logfile, *statfile;
FILE *coordfile, *restartcoordfile, *restartvelfile;
FILE *totdurationfile, *autodiffusionfile;
FILE *corrfuncfile;

//FILE *output, *cond_in, *cond_in_vel, *logfile, *param, *file_durata_totale;
//FILE *inputr, *inputv, *stat, *output_in_cond, *output_in_cond_vel;


vec rij, sumv, vi, rni;


float GaussianNoise(float mu, float var);
void set_initial_conditions(void);
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
void writePointCloudToVTK(const char *filename, const vec *points, int numPoints);
int readPointCloudFromVTK(const char *filename, vec *points, int numPoints);
void writePointCloudToVTKBinary(const char *filename, const vec *points, int numPoints);

