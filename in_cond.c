#include "defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

REAL GaussianNoise(REAL mu, REAL var){
  double u1, u2;
  u1 = rand() * (1.0 / RAND_MAX);
  u2 = rand() * (1.0 / RAND_MAX);
  double z0;
  //double z0, z1;
  z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
  //z1 = sqrt(-2.0 * log(u1)) * sin(2 * M_PI * u2);
  return z0 * var + mu;
}

void write_params(params p) {
  printf("twodim:         %d\n",p.twodim);
  printf("npartx:         %d\n",p.npartx);
  printf("nparty:         %d\n",p.nparty);
  printf("nlayers:        %d\n",p.nlayers);
  printf("npart:          %d\n",p.npart);
  printf("write_jump:     %d\n",p.write_jump);
  printf("timesteps:      %d\n",p.timesteps);
  printf("dt:             %lf\n",p.dt);
  printf("eps:            %lf\n",p.eps);
  printf("sigma:          %lf\n",p.sigma);
  printf("mu:             %lf\n",p.mu);
  printf("var:            %lf\n",p.var);
  printf("m:              %lf\n",p.m);
  printf("a_lattice:      %lf\n",p.a_lattice);
  printf("pot_trunc_perc: %lf\n",p.pot_trunc_perc);
  printf("newc:           %d\n",p.newc);
  printf("reproducible:   %d\n",p.reproducible);
  printf("r_max :         %lf\n",p.r_max);
  printf("BOXL:           %lf\n",p.BOXL);
  printf("shift:          %lf\n",p.shift);
  printf("dtsquare:       %lf\n",p.dtsquare);
  printf("dtdouble:       %lf\n",p.dtdouble);
  printf("r_max_squared:  %lf\n",p.r_max_squared);

  fprintf(p.logfile, "twodim:         %d\n",p.twodim);
  fprintf(p.logfile, "npartx:         %d\n",p.npartx);
  fprintf(p.logfile, "nparty:         %d\n",p.nparty);
  fprintf(p.logfile, "nlayers:        %d\n",p.nlayers);
  fprintf(p.logfile, "npart:          %d\n",p.npart);
  fprintf(p.logfile, "write_jump:     %d\n",p.write_jump);
  fprintf(p.logfile, "timesteps:      %d\n",p.timesteps);
  fprintf(p.logfile, "dt:             %lf\n",p.dt);
  fprintf(p.logfile, "eps:            %lf\n",p.eps);
  fprintf(p.logfile, "sigma:          %lf\n",p.sigma);
  fprintf(p.logfile, "mu:             %lf\n",p.mu);
  fprintf(p.logfile, "var:            %lf\n",p.var);
  fprintf(p.logfile, "m:              %lf\n",p.m);
  fprintf(p.logfile, "a_lattice:      %lf\n",p.a_lattice);
  fprintf(p.logfile, "pot_trunc_perc: %lf\n",p.pot_trunc_perc);
  fprintf(p.logfile, "newc:           %d\n",p.newc);
  fprintf(p.logfile, "reproducible:   %d\n",p.reproducible);
  fprintf(p.logfile, "r_max :         %lf\n",p.r_max);
  fprintf(p.logfile, "BOXL:           %lf\n",p.BOXL);
  fprintf(p.logfile, "shift:          %lf\n",p.shift);
  fprintf(p.logfile, "dtsquare:       %lf\n",p.dtsquare);
  fprintf(p.logfile, "dtdouble:       %lf\n",p.dtdouble);
  fprintf(p.logfile, "r_max_squared:  %lf\n",p.r_max_squared);
}

void set_initial_conditions(vec *r, vec *v, params p) {
  if (p.twodim==0) {
    if (p.npart != 4*p.nlayers*p.npartx*p.nparty) {
      printf("** PROBLEMA **\np.npart non congruente!\n");
      fprintf(p.logfile,"** PROBLEMA **\np.npart non congruente!\n");
    }
  } else {
    if (p.npart != 2*p.nlayers*p.npartx*p.nparty) {
      printf("** PROBLEMA **\np.npart non congruente!\n");
      fprintf(p.logfile,"** PROBLEMA **\np.npart non congruente!\n");
    }
  }
  REAL semi_lattice_len_x =  p.npartx*p.a_lattice/2.; // (0.5 + max(p.npartx-1,p.nparty-1,p.nlayers-1) * p.a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  REAL semi_lattice_len_y =  p.nparty*p.a_lattice/2.;
  REAL semi_lattice_len_z = p.nlayers*p.a_lattice/2.;
  
  int n = 0;
  if (p.twodim==0) {
    /**** To create FCC 100 lattice*********/
    for(int i = 0; i < p.npartx; i++){
      for(int j = 0; j < p.nparty; j++){
        for(int k = 0; k < p.nlayers; k++) { // Number of atoms in the Z direction
          r[n].x = i * p.a_lattice - semi_lattice_len_x;
          r[n].y = j * p.a_lattice - semi_lattice_len_y;
          r[n].z = k * p.a_lattice - semi_lattice_len_z;
          n++;

          r[n].x = i * p.a_lattice - semi_lattice_len_x;
          r[n].y = 0.5 * p.a_lattice +  j * p.a_lattice - semi_lattice_len_y;
          r[n].z = 0.5 * p.a_lattice + k * p.a_lattice - semi_lattice_len_z;
          n++;

          r[n].x = 0.5 * p.a_lattice + i * p.a_lattice - semi_lattice_len_x;
          r[n].y = j * p.a_lattice - semi_lattice_len_y;
          r[n].z = 0.5 * p.a_lattice + k *p.a_lattice - semi_lattice_len_z;
          n++;

          r[n].x = 0.5 * p.a_lattice + i*p.a_lattice - semi_lattice_len_x;
          r[n].y = 0.5 * p.a_lattice + j * p.a_lattice - semi_lattice_len_y;
          r[n].z = k * p.a_lattice - semi_lattice_len_z;
          n++;
        }
      }
    }
  } else {
    /**** To create 2d FCC lattice*********/
    for(int i = 0; i < p.npartx; i++){
      for(int j = 0; j < p.nparty; j++){
        r[n].x = (i+0.5) * p.a_lattice - semi_lattice_len_x;
        r[n].y = (j+0.5) * p.a_lattice - semi_lattice_len_y;
        r[n].z = 0.;
        n++;
        r[n].x = (i+1) * p.a_lattice - semi_lattice_len_x;
        r[n].y = (j+1) * p.a_lattice - semi_lattice_len_y;
        r[n].z = 0.;
        n++;
      }
    }

  }

  if (p.reproducible) {
    srand(0);
  } else {
    srand(time(NULL));
    }
  n = 0;
  for (int i = 0; i < p.npart; i++) { // p.npart / 2 perchè per annullare il momento e il momento angolare totale metto le velocità opposte a coppie
    v[n].x = GaussianNoise(p.mu, p.var);
    v[n].y = GaussianNoise(p.mu, p.var);
    v[n].z = GaussianNoise(p.mu, p.var);
    if (p.twodim==1) {v[n].z = 0.;}
    n++;
    v[n].x = -v[n-1].x; // giochetto velocità opposte a coppie
    v[n].y = -v[n-1].y;
    v[n].z = -v[n-1].z;
    if (p.twodim==1) {v[n].z = 0.;}
    n++;
  }
  
  // shuffle the array using Fisher-Yates algorithm
  for (int i = p.npart-1; i > 0; --i) {
    // Generate a random index j such that 0 <= j <= i
    int j = rand()%(i + 1);
    // Swap arr[i] with arr[j]
    vec tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
  }
}


void load_r_2(FILE *inputhere, vec *r, params p) {
  char garb[STRLEN];
  fscanf(inputhere,"%s\n\n",garb);
  for (int i = 0; i < p.npart; i++) {
    fscanf(inputhere,"atomX %lg %lg %lg\n", &(r[i].x), &(r[i].y), &(r[i].z));
  }
}


void write_stat(REAL t, vec sumv, REAL kenergy, REAL penergy, params p) {
  kenergy = kenergy * 0.5 * p.m;
  REAL reduced_temperature = 2. * kenergy / ((p.npart - 3.) * p.eps);
  fprintf(p.statfile,"%lg %lg %lg %lg %lg %lg %lg %lg\n", t*p.dt, sumv.x, sumv.y, sumv.z, kenergy, penergy, kenergy + penergy, reduced_temperature);
  printf("  tot mom: %lg %lg %lg\n  energy: %lg %lg %lg\n", sumv.x, sumv.y, sumv.z, kenergy, penergy, kenergy + penergy);
}


// Invio-numero particelle- invio (faccio capire che non si tratta di altre particelle,
// ma ho solo scalato lo step temporale)
void write_r(FILE *output_file, vec *r, params p) {
  fprintf(output_file, "%i\n\n", p.npart);
  for (int i = 0; i < p.npart; i++) {
    fprintf(output_file, "atomX %lg %lg %lg\n", r[i].x, r[i].y, r[i].z);
  }
}


void write_vi(FILE *output_file, vec vi) {
  fprintf(output_file, "%lg %lg %lg\n", vi.x, vi.y, vi.z);
}
