#include "defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

float GaussianNoise(float mu, float var){
  double u1, u2;
  u1 = rand() * (1.0 / RAND_MAX);
  u2 = rand() * (1.0 / RAND_MAX);
  double z0, z1;
  z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(2 * M_PI * u2);
  return z0 * var + mu;
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
  float semi_lattice_len_x =  p.npartx*p.a_lattice/2.; // (0.5 + max(p.npartx-1,p.nparty-1,p.nlayers-1) * p.a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  float semi_lattice_len_y =  p.nparty*p.a_lattice/2.;
  float semi_lattice_len_z = p.nlayers*p.a_lattice/2.;
  
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
  //
  //// shuffle the array using Fisher-Yates algorithm
  //for (int i = p.npart-1; i > 0; --i) {
  //  // Generate a random index j such that 0 <= j <= i
  //  int j = rand()%(i + 1);
  //  // Swap arr[i] with arr[j]
  //  vec tmp = v[i];
  //  v[i] = v[j];
  //  v[j] = tmp;
  //}
}


void load_r_2(FILE *inputhere, vec *r, params p) {
  char garb[STRLEN];
  fscanf(inputhere,"%s\n\n",garb);
  for (int i = 0; i < p.npart; i++) {
    fscanf(inputhere,"atomX %g %g %g\n", &(r[i].x), &(r[i].y), &(r[i].z));
  }
}


void write_stat(float t, vec sumv, float kenergy, float penergy, params p) {
  kenergy = kenergy * 0.5 * p.m;
  float reduced_temperature = 2. * kenergy / ((p.npart - 3.) * p.eps);
  fprintf(p.statfile,"%g %g %g %g %g %g %g %g\n", t*p.dt, sumv.x, sumv.y, sumv.z, kenergy, penergy, kenergy + penergy, reduced_temperature);
  printf("tot mom: %g %g %g\nenergy: %g %g %g\n", sumv.x, sumv.y, sumv.z, kenergy, penergy, kenergy + penergy);
}


// Invio-numero particelle- invio (faccio capire che non si tratta di altre particelle,
// ma ho solo scalato lo step temporale)
void write_r(FILE *output_file, vec *r, params p) {
  fprintf(output_file, "%i\n\n", p.npart);
  for (int i = 0; i < p.npart; i++) {
    fprintf(output_file, "atomX %g %g %g\n", r[i].x, r[i].y, r[i].z);
  }
}


void write_vi(FILE *output_file, vec vi) {
  fprintf(output_file, "%g %g %g\n", vi.x, vi.y, vi.z);
}
