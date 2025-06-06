#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

#define CHECK_FILE(ptr, name) \
  if (!(ptr)) { perror("fopen " #name); exit(1); }

int main(int argc, char *argv[]) {
  // CHECK IF COMMAND LINE ARGUMENTS ARE PROVIDED
  if (argc < 2) {fprintf(stderr, "Usage: %s <input_file_path>\n", argv[0]); return -1;}
  // WRITE INPUT AND OUTPUT PATHS
  snprintf(paramfilepath,        STRLEN, "%s/param.in",            argv[1]);
  snprintf(logfilepath,          STRLEN, "%s/run.log",             argv[1]);
  snprintf(statfilepath,         STRLEN, "%s/stat.dat",            argv[1]);
  snprintf(coordfilepath,        STRLEN, "%s/verlet_periodic.xyz", argv[1]);
  snprintf(restartcoordfilepath, STRLEN, "%s/in_cond.xyz",         argv[1]);
  snprintf(restartvelfilepath,   STRLEN, "%s/in_cond_vel.dat",     argv[1]);
  snprintf(restartvelfilepath,   STRLEN, "%s/in_cond_vel.dat",     argv[1]);
  snprintf(totdurationfilepath,  STRLEN, "%s/tot_duration.dat",     argv[1]);
  paramfile        = fopen(paramfilepath, "r");
  logfile          = fopen(logfilepath, "w");
  statfile         = fopen(statfilepath, "w");
  coordfile        = fopen(coordfilepath, "w");
  restartcoordfile = fopen(restartcoordfilepath, "w");
  restartvelfile   = fopen(restartvelfilepath, "w");
  totdurationfile  = fopen(totdurationfilepath, "w");
  CHECK_FILE(paramfile,        paramfilepath);
  CHECK_FILE(logfile,          logfilepath);
  CHECK_FILE(statfile,         statfilepath);
  CHECK_FILE(coordfile,        coordfilepath);
  CHECK_FILE(restartcoordfile, restartcoordfilepath);
  CHECK_FILE(restartvelfile,   restartvelfilepath);
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

  if (newc==0) {printf("Restart simulation feature not available. Stoppingi\n"); return -1;}
  printf("r_max=%g    BOXL=%g    red. dens=%g\n",r_max,BOXL,reduced_density);
  nrun = 0;
  last_durata_totale = 0.;
  printf("Initialize FCC lattice and random velocities.\n");
  fprintf(logfile,"start new simulation:\n\nPARAM:\n");
  fprintf(logfile,"npartx=%i\nnparty=%i\n",npartx,nparty);
  fprintf(logfile,"nlayers=%i\nnpart=%i\n",nlayers,npart);
  fprintf(logfile,"write_jump=%i\ntimesteps=%i\n",write_jump,timesteps);
  fprintf(logfile,"dt=%g\neps=%g\nsigma=%g\n",dt,eps,sigma);
  fprintf(logfile,"mu=%g\nvar=%g\nm=%g\na_lattice=%g\n",mu,var,m,a_lattice);
  fprintf(logfile,"pot_trunc_perc=%g\nnew_in_cond=%i\n\n\n",pot_trunc_perc,newc);
  fprintf(logfile,"Initialize FCC lattice and random velocities\n\n");
  fprintf(logfile, "r_max=%g    BOXL=%g    red. dens=%g\n",r_max,BOXL,reduced_density);
  fcc();
    
    
  // OPEN OUTPUT FILES AND WRITE HEADERS
  fprintf(statfile,"# t(0)   sumvx(1)   sumvy(2)   sumvz(3)    kenergy(4)   penergy(5)   energy(6)   red_temp(7)\n");
    
  // DECLARE VARIABLES
  vec r[npart], ro[npart], a[npart];
  
  // LOAD INITIAL CONDITIONS
  load_r(r);
  write_r(coordfile, r);
  printf("Integration started:\n\n");
  fprintf(logfile, "Integration started: wait!\n\n");
    for (int i = 0; i < npart; i++) {
        a[i].x = 0.0;
        a[i].y = 0.0;
        a[i].z = 0.0;
    }
  compute_forces(r, a);
  eulero(r, ro, a);
  write_r(coordfile, r);


  // CHECK PARTICLES INSIDE BOX
  int out_count = 0;
  for (int i = 0; i < npart; i++) {
    if (r[i].x < -BOXL/2. || r[i].x > BOXL/2.) out_count += 1;
    if (r[i].y < -BOXL/2. || r[i].y > BOXL/2.) out_count += 1;
    if (r[i].z < -BOXL/2. || r[i].z > BOXL/2.) out_count += 1;
  }
  if (out_count > 0) printf("initial conditions are not inside periodic box.\n");


  // VERLET INTEGRATION
  t = 2;
  for (int tt = 2; tt < (timesteps-1) / write_jump; tt++) {
    for (int k = 0; k < write_jump; k++) {
      compute_forces(r, a);
      verlet_periodic(r, ro, a);
      t++;
    }
    compute_forces_stat(r, a);
    verlet_periodic_write(r, ro, a);
    t++;
  }
  
  // LAST TIME STEP
  compute_forces_stat(r, a);
  verlet_periodic_last(r, ro, a); t++;
  write_durata_totale();
  
  // CLOSE FILES
  fclose(coordfile);
  fclose(restartcoordfile);
  fclose(restartvelfile);
  fclose(statfile);
  printf("Done!\n\n");
  fprintf(logfile, "Done!\n\n\n\n\n\n\n\n");
  fclose(logfile);
  fclose(totdurationfile);
}

/*
NOTE:
1) Ho notato che ponendo pot_tronc_perc = 1-e-3, ovvero r_max \simeq 3*sigma, si introduce un piccolo errore sulla funzione
  dell'energia cinetica e potenziale in funzione del tempo (e probabilmente anche su tutto il resto, tipo le coordinate
  delle particelle). Per indagare ulteriormente ciò andrebbero analizzati gli errori relativi tra il potenziale non troncato
  e quello troncato con vari valori di pot_tronc_perc... Fino ad arrivare a 1.e-7, nel quale pare che il profili delle energie
  siano "uguali".

2) Shiftando solo il potenziale con una costante ho introdotto una discostinuità nella forza (diventa nulla all'improvviso);
   Potrei correggere questa cosa inserendo un termine lineare piccolo che fa diventare la forza nulla con continuità quando
  rij = r_max, ma introducendo nella forza un termine lineare in rij poi dovrei calcolare la radice quadrata e sarebbe molto
  più inefficiente.

*/
