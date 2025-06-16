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


void set_initial_conditions(params p) {
  if (p.npart != 4*p.nlayers*p.npartx*p.nparty) {
    printf("** PROBLEMA **\np.npart non congruente!\n");
    fprintf(p.logfile,"** PROBLEMA **\np.npart non congruente!\n");
  }
  vec r, v0; // giochetto velocità opposte a coppie (mom angolare e mom lineare circa nulli)
  float vxm = 0., vym = 0., vzm = 0.;
  char AtomName[] = "atomX";
  float semi_lattice_len_x =  p.npartx*p.a_lattice/2.; // (0.5 + max(p.npartx-1,p.nparty-1,p.nlayers-1) * p.a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  float semi_lattice_len_y =  p.nparty*p.a_lattice/2.;
  float semi_lattice_len_z = p.nlayers*p.a_lattice/2.;
  
  FILE *cond_in = fopen(p.restartcoordfilepath,"w");
  if (!cond_in) {
    perror("fopen cond_in");
    exit(1);
  }
  FILE *cond_in_vel = fopen(p.restartvelfilepath,"w");
  if (!cond_in_vel) {
    perror("fopen cond_in_vel");
    exit(1);
  }

  fprintf(cond_in,"%i\n\n",p.npart);

  /**** To create FCC 100 lattice*********/
  for(int i = 0; i < p.npartx; i++){            // Number of Atoms in the X direction

    for(int j = 0; j < p.nparty; j++){        // Number of Atoms in the Y direction

      for(int k = 0; k < p.nlayers; k++) {       // Number of layers in the Z direction

        r.x = i * p.a_lattice - semi_lattice_len_x;   r.y = j * p.a_lattice - semi_lattice_len_y;   r.z = k * p.a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);

        r.x = i * p.a_lattice - semi_lattice_len_x;   r.y = 0.5 * p.a_lattice +  j * p.a_lattice - semi_lattice_len_y;    r.z = 0.5 * p.a_lattice + k * p.a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);

        r.x = 0.5 * p.a_lattice + i * p.a_lattice - semi_lattice_len_x;  r.y = j * p.a_lattice - semi_lattice_len_y;   r.z = 0.5 * p.a_lattice + k *p.a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);

        r.x = 0.5 * p.a_lattice + i*p.a_lattice - semi_lattice_len_x;  r.y = 0.5 * p.a_lattice + j * p.a_lattice - semi_lattice_len_y;   r.z = k * p.a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);
      }

     }

  }

  if (p.reproducible) {
    srand(0);
  } else {
    srand(time(NULL));
    }
  for (int i = 0; i < p.npart / 2; i++) { // p.npart / 2 perchè per annullare il momento e il momento angolare totale metto le velocità opposte a coppie
    v0.x = GaussianNoise(p.mu, p.var);
    v0.y = GaussianNoise(p.mu, p.var);
    v0.z = GaussianNoise(p.mu, p.var);
    fprintf(cond_in_vel,"%g %g %g\n", v0.x, v0.y, v0.z); // giochetto velocità opposte a coppie
    fprintf(cond_in_vel,"%g %g %g\n",-v0.x,-v0.y,-v0.z);
  }
  fclose(cond_in);
  fclose(cond_in_vel);
}


void load_v(vec *r, params p) {
  FILE *inputr = fopen(p.restartvelfilepath,"r");
  for (int i = 0; i < p.npart; i++) {
    fscanf(inputr,"%g %g %g\n", &(r[i].x), &(r[i].y), &(r[i].z));
  }
  fclose(inputr);
}


void load_r(vec *r, params p) {
  FILE *inputr = fopen(p.restartcoordfilepath,"r");
  char garb[STRLEN];
  fscanf(inputr,"%s\n\n",garb);
  for (int i = 0; i < p.npart; i++) {
    fscanf(inputr,"atomX %g %g %g\n", &(r[i].x), &(r[i].y), &(r[i].z));
        //printf("%g %g %g",r[i].x,r[i].y,r[i].z);
  }
  fclose(inputr);
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
