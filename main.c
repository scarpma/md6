#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include <time.h>

#define CHECK_FILE(ptr, name) \
  if (!(ptr)) { perror("fopen " #name); exit(1); }

int main(int argc, char *argv[]) {
  params p;
  int t;
  REAL penergy;
  printf("Ciao dal corso, (non) sono Filippo. Chi sono?\n");
  printf("Hola senor Martino 💃💃💃💃💃💃💃!\n");
  // CHECK IF COMMAND LINE ARGUMENTS ARE PROVIDED
  if (argc < 2) {fprintf(stderr, "Usage: %s <input_file_path>\n", argv[0]); return -1;}
  // WRITE INPUT AND OUTPUT PATHS
  snprintf(p.paramfilepath,        STRLEN, "%s/param.in",            argv[1]);
  snprintf(p.logfilepath,          STRLEN, "%s/run.log",             argv[1]);
  snprintf(p.statfilepath,         STRLEN, "%s/stat.dat",            argv[1]);
  snprintf(p.coordfilepath,        STRLEN, "%s/verlet_periodic.xyz", argv[1]);
  snprintf(p.restartcoordfilepath, STRLEN, "%s/in_cond.xyz",         argv[1]);
  snprintf(p.restartvelfilepath,   STRLEN, "%s/in_cond_vel.dat",     argv[1]);
  snprintf(p.restartvelfilepath,   STRLEN, "%s/in_cond_vel.dat",     argv[1]);
  p.paramfile        = fopen(p.paramfilepath, "r");
  p.logfile          = fopen(p.logfilepath, "w");
  p.statfile         = fopen(p.statfilepath, "w");
  p.coordfile        = fopen(p.coordfilepath, "w");
  p.restartcoordfile = fopen(p.restartcoordfilepath, "w");
  p.restartvelfile   = fopen(p.restartvelfilepath, "w");
  CHECK_FILE(p.paramfile,        p.paramfilepath);
  CHECK_FILE(p.logfile,          p.logfilepath);
  CHECK_FILE(p.statfile,         p.statfilepath);
  CHECK_FILE(p.coordfile,        p.coordfilepath);
  CHECK_FILE(p.restartcoordfile, p.restartcoordfilepath);
  CHECK_FILE(p.restartvelfile,   p.restartvelfilepath);
  char vtkfilename[STRLEN];
  fprintf(p.statfile,"# t(0)   sumvx(1)   sumvy(2)   sumvz(3)    kenergy(4)   penergy(5)   energy(6)   red_temp(7)\n");

  // READ PARAMETERS
  fscanf(p.paramfile,"twodim=%d\nnpartx=%d\nnparty=%d\nnlayers=%d\nnpart=%d\nwrite_jump=%d\ntimesteps=%d\ndt=%lg\neps=%lg\nsigma=%lg\nmu=%lg\nvar=%lg\nm=%lg\na_lattice=%lg\npot_trunc_perc=%lg\nnew_in_cond=%d\nreproducible=%d",&p.twodim,&p.npartx,&p.nparty,&p.nlayers,&p.npart,&p.write_jump,&p.timesteps,&p.dt,&p.eps,&p.sigma,&p.mu,&p.var,&p.m,&p.a_lattice,&p.pot_trunc_perc,&p.newc,&p.reproducible);
  fclose(p.paramfile);
  // INITIALIZE VARIABLES
  p.dtdouble = 2.*p.dt;
  p.dtsquare = pow(p.dt, 2.);
  p.r_max = p.sigma * pow((1+sqrt( 1-16*p.pot_trunc_perc ))/( 2*p.pot_trunc_perc ), 1./6.);  /*r_max IS COMPUTED BASED ON pot_trunc_perc, i.e. when POTENTIAL REACHES pot_trunc_perc OF ITS MAX VALUE*/
  p.r_max_squared = pow(p.r_max, 2.);
  p.shift = potenergy(p.r_max, p.eps, p.sigma); // POTENTIAL SHIFT
  if (p.twodim==0) {
    p.BOXL = ( 0.5 + p.nlayers ) * p.a_lattice; // (0.5 + max(npartx,nparty,nlayers) * a_lattice IS THE LATTICE LENGHT IN EACH DIRECTION)
  } else {
    p.BOXL = p.npartx * p.a_lattice;
  }
  p.reduced_density = p.npart * p.sigma / pow(p.BOXL, 3.);

  if (p.newc==0) {printf("Restart simulation feature not available. Stoppingi\n"); return -1;}
  write_params(p);
    
  // DECLARE VARIABLES
  vec r[p.npart], ro[p.npart], v[p.npart], a[p.npart];
  for (int i=0; i<p.npart; i++) {
    r[i].x  = 0.;
    r[i].y  = 0.;
    r[i].z  = 0.;
    ro[i].x = 0.;
    ro[i].y = 0.;
    ro[i].z = 0.;
    v[i].x  = 0.;
    v[i].y  = 0.;
    v[i].z  = 0.;
    a[i].x  = 0.;
    a[i].y  = 0.;
    a[i].z  = 0.;
  }
  set_initial_conditions(r, v, p);
    
  
  printf("Integration started:\n\n");
  fprintf(p.logfile, "Integration started: wait!\n\n");

  printf("Ciao dal corso");
  printf("from the past");

  // CHECK PARTICLES INSIDE BOX
  int out_count = 0;
  for (int i = 0; i < p.npart; i++) {
    if (r[i].x < -p.BOXL/2. || r[i].x > p.BOXL/2.) out_count += 1;
    if (r[i].y < -p.BOXL/2. || r[i].y > p.BOXL/2.) out_count += 1;
    if (r[i].z < -p.BOXL/2. || r[i].z > p.BOXL/2.) out_count += 1;
  }
  if (out_count > 0) printf("initial conditions are not inside periodic box.\n");


  // VERLET INTEGRATION
  clock_t start, end;
  t = 0;
  while (t < p.timesteps) {
    if (t%p.write_jump==0) {
      start = clock();
      snprintf(vtkfilename, STRLEN, "%s/particles_%08d.vtk", argv[1], t);
      //writePointCloudToVTKBinary(vtkfilename, r, p.npart); // does not work
      writePointCloudToVTKBinaryBigEndian(vtkfilename, r, v, p.npart);
      //writePointCloudToVTK(vtkfilename, r, v, p.npart);
    }


    // =========================== //
    // POTENTIAL ENERGY AND FORCES //
    // =========================== //
    if (t%p.write_jump==0) {
      penergy = compute_forces_stat(r, a, p);
    } else {
      compute_forces(r, a, p);
    }

    if (p.twodim==1) {for (int i = 0; i < p.npart; i++) {r[i].z=0.;}}


    // =========================== //
    // KINETIC ENERGY AND MOMENTUM //
    // =========================== //
    
    if (t%p.write_jump==0) {
      if (t==0) {
        compute_kenergy_momentum2(t, r, v, penergy, p);
      } else {
        compute_kenergy_momentum(t, r, ro, v, penergy, p);
      }
    }


    // =========================== //
    // NEWTON MECHANICS TIMESTEP   //
    // =========================== //
    if (t==0) {
      eulero(r, ro, v, a, p);
    } else {
      verlet_periodic(r, ro, a, p);
    }


    if (t%p.write_jump==0) {
      end = clock();
      double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
      printf("Iteration: %d, %.3f msec\n", t, time_taken*1000);
    }


    t++;
  }
  
  // CLOSE FILES
  fclose(p.coordfile);
  fclose(p.restartcoordfile);
  fclose(p.restartvelfile);
  fclose(p.statfile);
  printf("Done!\n\n");
  fprintf(p.logfile, "Done!\n\n\n\n\n\n\n\n");
  fclose(p.logfile);
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
