#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// useful for debbugging
#define LOG_LOCATION() printf("File: %s, Line: %d\n", __FILE__, __LINE__)
#define STRLEN 256

struct vec3d {float x, y, z;};
typedef struct vec3d vec;

struct paramsMD {
  // read params
  int npartx;
  int nparty;
  int nlayers;
  int npart;
  int write_jump;
  int timesteps;
  float dt;
  float eps;
  float sigma;
  float mu;
  float var;
  float m;
  float a_lattice;
  float pot_trunc_perc;
  int newc;
  int reproducible;

  // computed params
  float r_max, BOXL, shift;
  float reduced_density, red_temp;
  // handy variables used often
  float dtsquare, dtdouble, r_max_squared;
  
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

float GaussianNoise(float mu, float var);
void set_initial_conditions(params p);
void load_r(vec *r, params p);
void load_r_2(FILE *inputhere, vec *r, params p);
void write_stat(float t, vec sumv, float kenergy, float penergy, params p);
void write_r(FILE *output_file, vec *r, params p);
void write_vi(FILE *output_file, vec vi);

float potenergy(float r, float eps, float sigma);
float simforce(float r, float eps, float sigma);
void verlet_periodic(vec *r, vec *ro, vec *a, params p);
void compute_kenergy_momentum(float t, vec *r, vec *ro, vec *a, float penergy, params p);
void compute_forces(vec *r, vec *a, params p);
float compute_forces_stat(vec *r, vec *a, params p);
void eulero(vec *r, vec *ro, vec *a, params p);
void writePointCloudToVTK(const char *filename, const vec *points, int numPoints);
void writePointCloudToVTKBinary(const char *filename, const vec *points, int numPoints);

