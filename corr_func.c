#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

// Define variables once and for all
int reproducible;
int t, npartx, nparty, nlayers, npart, write_jump, timesteps, newc, nrun, garb;
float dt, dtsquare, dtdouble, eps, sigma, mu, var, m, a_lattice, last_durata_totale;
float r_max, pot_trunc_perc, BOXL, r_max_squared, shift, penergy, rij2;
float kenergy, simforceij, reduced_temperature, reduced_density;

char paramfilepath[STRLEN], logfilepath[STRLEN];
char statfilepath[STRLEN];
char coordfilepath[STRLEN], restartcoordfilepath[STRLEN];
char restartvelfilepath[STRLEN], autodiffusionfilepath[STRLEN];
char totdurationfilepath[STRLEN], corrfuncfilepath[STRLEN];

FILE *paramfile, *logfile, *statfile;
FILE *coordfile, *restartcoordfile, *restartvelfile;
FILE *totdurationfile, *autodiffusionfile;
FILE *corrfuncfile;

vec rij, sumv, vi, rni;

#define CHECK_FILE(ptr, name) \
  if (!(ptr)) { perror("fopen " #name); exit(1); }

int main(int argc, char *argv[]) {
  // CHECK IF COMMAND LINE ARGUMENTS ARE PROVIDED
  if (argc < 2) {fprintf(stderr, "Usage: %s <input_file_path>\n", argv[0]); return -1;}
  // APRO FILES E LI INIZIALIZZO
  snprintf(paramfilepath,         STRLEN, "%s/param.in",            argv[1]);
  snprintf(coordfilepath,         STRLEN, "%s/verlet_periodic.xyz", argv[1]);
  snprintf(corrfuncfilepath,      STRLEN, "%s/corr_func.dat",   argv[1]);
  snprintf(totdurationfilepath,   STRLEN, "%s/tot_duration.dat",    argv[1]);
  paramfile        = fopen(paramfilepath,            "r");
  coordfile        = fopen(coordfilepath,            "r");
  corrfuncfile     = fopen(corrfuncfilepath,         "w");
  totdurationfile   = fopen(totdurationfilepath, "r");
  CHECK_FILE(paramfile,        paramfilepath);
  CHECK_FILE(coordfile,        coordfilepath);
  CHECK_FILE(corrfuncfile,     corrfuncfilepath);
  CHECK_FILE(totdurationfile,  totdurationfilepath);

  // READ PARAMETERS
  fscanf(paramfile,"npartx=%i\nnparty=%i\nnlayers=%i\nnpart=%i\nwrite_jump=%i\ntimesteps=%i\ndt=%g\neps=%g\nsigma=%g\nmu=%g\nvar=%g\nm=%g\na_lattice=%g\npot_trunc_perc=%g\nnew_in_cond=%i",&npartx,&nparty,&nlayers,&npart,&write_jump,&timesteps,&dt,&eps,&sigma,&mu,&var,&m,&a_lattice,&pot_trunc_perc,&newc);
  fclose(paramfile);

  // INITIALIZE VARIABLES
  dtdouble = 2.*dt;
  dtsquare = pow(dt, 2.);
  r_max = sigma * pow((1+sqrt( 1-16*pot_trunc_perc ))/( 2*pot_trunc_perc ), 1./6.);  /*r_max IS COMPUTED BASED ON pot_trunc_perc, i.e. when POTENTIAL REACHES pot_trunc_perc OF ITS MAX VALUE*/
  r_max_squared = pow(r_max, 2.);
  shift = potenergy(r_max, eps, sigma); // POTENTIAL SHIFT
  BOXL = ( 0.5 + nlayers ) * a_lattice; // (0.5 + max(npartx,nparty,nlayers) * a_lattice IS THE LATTICE LENGHT IN EACH DIRECTION)
  reduced_density = npart * sigma / pow(BOXL, 3.);

  // CARICO DURATA TOTALE
  fscanf(totdurationfile,"%g\n%i\n",&last_durata_totale,&nrun);
  
  // INIZIALIZZO VARIABILI AUTODIFFUSION
  int ntime = timesteps / write_jump;
  float dx = sigma / 30.;
  int n_bin = (int)(r_max / dx);
  float xv[n_bin], fv[n_bin];
  for (int a = 0; a < n_bin; a++) xv[a] = dx*a + 0.5*dx;
  for (int a = 0; a < n_bin; a++) fv[a] = 0;
  vec r[npart], rij;
  float rij2;

  // INIZIALIZZO FILE OUTPUT
  fprintf(corrfuncfile, "#    x(1)      f(x)(2)\n");

  // CARICO r
  for (int k = 0; k < ntime; k++) {
  load_r_2(coordfile, r);
    
    // CICLO SU OGNI PARTICELLA PRIMARIA j
    for (int j = 0; j < npart; j++) {
      
      // CICLO SU OGNI PARTICELLA SECONDARIA i
      for (int i = 0; i < npart; i++) {
        rij.x = (r[j].x - r[i].x) - BOXL*round((r[j].x-r[i].x)/BOXL);
        rij.y = (r[j].y - r[i].y) - BOXL*round((r[j].y-r[i].y)/BOXL);
        rij.z = (r[j].z - r[i].z) - BOXL*round((r[j].z-r[i].z)/BOXL);
        rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
        if (rij2 < r_max_squared && rij2 != 0.) {
          rij2 = pow(rij2, 0.5);
          fv[(int)floor(rij2/dx)]++;
        }
      }
    }
  }
  for (int a = 0; a < n_bin; a++) {
    fv[a] /= (npart*ntime);
    fprintf(corrfuncfile, "%g %g %g\n",xv[a],fv[a],pow(xv[a],2.)*4.*3.1415*dx);
  }
  fclose(corrfuncfile);
}

