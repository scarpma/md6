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
  snprintf(autodiffusionfilepath, STRLEN, "%s/autodiffusion.dat",   argv[1]);
  snprintf(totdurationfilepath,   STRLEN, "%s/tot_duration.dat",    argv[1]);
  paramfile         = fopen(paramfilepath,            "r");
  coordfile         = fopen(coordfilepath,            "r");
  autodiffusionfile = fopen(autodiffusionfilepath, "w");
  totdurationfile   = fopen(totdurationfilepath, "r");
  CHECK_FILE(paramfile,        paramfilepath);
  CHECK_FILE(coordfile,        coordfilepath);
  CHECK_FILE(autodiffusionfile,autodiffusionfilepath);
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
  vec r[npart], rin[npart], disp2[npart], var;
  float time;

  // INIZIALIZZO FILE OUTPUT
  fprintf(autodiffusionfile, "#    t(1)      varx(2)      vary(3)      varz(3)\n");
  
  // CICLO SUI TEMPI

  for (int j = 0; j < ntime; j++) { // ciclo su ogni istante di tempo (scritto nel file)

    // CARICO r e rin
    load_r_2(coordfile, r);
    // for (int i = 0; i < npart; i++) printf("%g %g %g\n",r[i].x,r[i].y,r[i].z);
    if (j == 0) {  
      for (int i = 0; i < npart; i++) {
        rin[i].x = r[i].x;
        rin[i].y = r[i].y;
        rin[i].z = r[i].z;
      }
    }


    // CICLO SULLE PARTICELLE AL TEMPO j (NO ! Non va bene così. Perchè comunque rx -rxin non sarà mai maggiore di ????)
    for (int i = 0; i < npart; i++) {
      disp2[i].x = pow(( r[i].x - rin[i].x ) - BOXL*round(( r[i].x - rin[i].x )/BOXL), 2.);
      disp2[i].y = pow(( r[i].y - rin[i].y ) - BOXL*round(( r[i].y - rin[i].y )/BOXL), 2.);
      disp2[i].z = pow(( r[i].z - rin[i].z ) - BOXL*round(( r[i].z - rin[i].z )/BOXL), 2.);
      //printf("%g %g %g\n",disp2[i].x,disp2[i].y,disp2[i].z);          
    }
    var.x = 0.; var.y = 0.; var.z = 0.;
    for (int i = 0; i < npart; i++) {
      var.x += disp2[i].x;
      var.y += disp2[i].y;
      var.z += disp2[i].z;
    }
    time = j * dt * write_jump + last_durata_totale - timesteps*dt;
    var.x = var.x/npart;
    var.y = var.y/npart;
    var.z = var.z/npart;
    fprintf(autodiffusionfile, "%g %g %g %g\n",time,var.x,var.y,var.z);

  }
}





