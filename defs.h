#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// useful for debbugging
#define LOG_LOCATION() printf("File: %s, Line: %d\n", __FILE__, __LINE__)
#define STRLEN 256
#define REAL double 

//FLOAT IS NOT SUPPORTED ANYMORE
//#define REAL float

struct vec3d {REAL x, y, z;};
typedef struct vec3d vec;

struct paramsMD {
  // read params
  int twodim;
  int npartx;
  int nparty;
  int nlayers;
  int npart;
  int write_jump;
  int timesteps;
  REAL dt;
  REAL eps;
  REAL sigma;
  REAL mu;
  REAL var;
  REAL m;
  REAL a_lattice;
  REAL pot_trunc_perc;
  int newc;
  int reproducible;

  // computed params
  REAL r_max, BOXL, shift;
  REAL reduced_density, red_temp;
  // handy variables used often
  REAL dtsquare, dtdouble, r_max_squared;
  
  // path of files
  char paramfilepath[STRLEN];
  char logfilepath[STRLEN];
  char statfilepath[STRLEN];
  char coordfilepath[STRLEN];
  char restartcoordfilepath[STRLEN];
  char restartvelfilepath[STRLEN];
  char autodiffusionfilepath[STRLEN];
  char corrfuncfilepath[STRLEN];
  
  // file pointers
  FILE *paramfile;
  FILE *logfile;
  FILE *statfile;
  FILE *coordfile;
  FILE *restartcoordfile;
  FILE *restartvelfile;
  FILE *autodiffusionfile;
  FILE *corrfuncfile;
};
typedef struct paramsMD params;

void write_params(params p);

REAL GaussianNoise(REAL mu, REAL var);
void set_initial_conditions(vec *r, vec *v, params p);
void load_r(vec *r, params p);
void load_v(vec *v, params p);
void load_r_2(FILE *inputhere, vec *r, params p);
void write_stat(REAL t, vec sumv, REAL kenergy, REAL penergy, params p);
void write_r(FILE *output_file, vec *r, params p);
void write_vi(FILE *output_file, vec vi);

REAL potenergy(REAL r, REAL eps, REAL sigma);
REAL simforce(REAL r, REAL eps, REAL sigma);
void verlet_periodic(vec *r, vec *ro, vec *a, params p);
void compute_kenergy_momentum(REAL t, vec *r, vec *ro, vec *v, REAL penergy, params p);
void compute_kenergy_momentum2(REAL t, vec *r, vec *v, REAL penergy, params p);
void compute_forces(vec *r, vec *a, params p);
REAL compute_forces_stat(vec *r, vec *a, params p);
void eulero(vec *r, vec *ro, vec *v, vec *a, params p);
void writePointCloudToVTK(const char *filename, const vec *points, const vec *velocity, int numPoints);
void writePointCloudToVTKBinary(const char *filename, const vec *points, int numPoints);
void writePointCloudToVTKBinaryBigEndian(const char *filename, const vec *points, const vec *velocity, int numPoints);

