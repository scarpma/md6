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


void set_initial_conditions(vec *r, vec *v0) {
  if (nlayers==1) {
    if (npartx*nparty!=npart) {
      printf("** PROBLEMA **\nnpart non congruente!\n");
      fprintf(logfile,"** PROBLEMA **\nnpart non congruente!\n");
      //exit(-1);
    }
  } else {
    if (npart != 4*nlayers*npartx*nparty) {
      printf("** PROBLEMA **\nnpart non congruente!\n");
      fprintf(logfile,"** PROBLEMA **\nnpart non congruente!\n");
      exit(-1);
    }
  }
  float vxm = 0., vym = 0., vzm = 0.;
  float semi_lattice_len_x =  npartx*a_lattice/2.; // (0.5 + max(npartx-1,nparty-1,nlayers-1) * a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  float semi_lattice_len_y =  nparty*a_lattice/2.;
  float semi_lattice_len_z = nlayers*a_lattice/2.;
  
  int n = 0;
  if (nlayers>1) {
    /**** To create FCC 100 lattice*********/
    for(int i = 0; i < npartx; i++){      // Number of Atoms in the X direction
      for(int j = 0; j < nparty; j++){    // Number of Atoms in the Y direction
        for(int k = 0; k < nlayers; k++) {// Number of layers in the Z direction

          r[n].x = i * a_lattice - semi_lattice_len_x;
          r[n].y = j * a_lattice - semi_lattice_len_y;
          r[n].z = k * a_lattice - semi_lattice_len_z;
          n++;

          r[n].x = i * a_lattice - semi_lattice_len_x;
          r[n].y = 0.5 * a_lattice +  j * a_lattice - semi_lattice_len_y;
          r[n].z = 0.5 * a_lattice + k * a_lattice - semi_lattice_len_z;
          n++;

          r[n].x = 0.5 * a_lattice + i * a_lattice - semi_lattice_len_x;
          r[n].y = j * a_lattice - semi_lattice_len_y;
          r[n].z = 0.5 * a_lattice + k *a_lattice - semi_lattice_len_z;
          n++;

          r[n].x = 0.5 * a_lattice + i*a_lattice - semi_lattice_len_x;
          r[n].y = 0.5 * a_lattice + j * a_lattice - semi_lattice_len_y;
          r[n].z = k * a_lattice - semi_lattice_len_z;
          n++;
        }
      }
    }
  } else {
    /**** To create 2d FCC lattice*********/
    int n = 0;
    for(int i = 0; i < npartx; i++){
      for(int j = 0; j < nparty; j++){
        r[n].x = (i+0.5) * a_lattice - semi_lattice_len_x;
        r[n].y = (j+0.5) * a_lattice - semi_lattice_len_y;
        r[n].z = 0.;
        n++;
        r[n].x = (i+1) * a_lattice - semi_lattice_len_x;
        r[n].y = (j+1) * a_lattice - semi_lattice_len_y;
        r[n].z = 0.;
        n++;
      }
    }

  }

  if (reproducible) {
    srand(0);
  } else {
    srand(time(NULL));
    }
  for (int i = 0; i < npart-1; i=i+2) {
    // npart / 2 perchè per annullare il momento
    // e il momento angolare totale metto le velocità
    // opposte a coppie
    v0[i].x = GaussianNoise(mu, var);
    v0[i].y = GaussianNoise(mu, var);
    v0[i].z = GaussianNoise(mu, var);
    v0[i+1].x = -v0[i].x;
    v0[i+1].y = -v0[i].y;
    v0[i+1].z = -v0[i].z;
  }
  //// shuffle the array using Fisher-Yates algorithm
  //for (int i = npart-1; i > 0; --i) {
  //  // Generate a random index j such that 0 <= j <= i
  //  int j = rand()%(i + 1);
  //  // Swap arr[i] with arr[j]
  //  vec tmp = v0[i];
  //  v0[i] = v0[j];
  //  v0[j] = tmp;
  //}
}


void load_r(vec *r) {
  FILE *inputr = fopen(restartcoordfilepath,"r");
  fscanf(inputr,"%i\n\n",&garb);
  for (int i = 0; i < npart; i++) {
    fscanf(inputr,"atomX %g %g %g\n", &(r[i].x), &(r[i].y), &(r[i].z));
        //printf("%g %g %g",r[i].x,r[i].y,r[i].z);
  }
  fclose(inputr);
}


void load_r_2(FILE *inputhere, vec *r) {
  fscanf(inputhere,"%i\n\n",&garb);
  for (int i = 0; i < npart; i++) {
    fscanf(inputhere,"atomX %g %g %g\n", &(r[i].x), &(r[i].y), &(r[i].z));
  }
}


void write_stat(void) {
  kenergy = kenergy * 0.5 * m;
  reduced_temperature = 2. * kenergy / ((npart - 3.) * eps);
  fprintf(statfile,"%g %g %g %g %g %g %g %g\n", t*dt + last_durata_totale, sumv.x, sumv.y, sumv.z, kenergy, penergy, kenergy + penergy, reduced_temperature);
}


// Invio-numero particelle- invio (faccio capire che non si tratta di altre particelle,
// ma ho solo scalato lo step temporale)
void write_r(FILE *output_file, vec *r) {
  fprintf(output_file, "%i\n\n", npart);
  for (int i = 0; i < npart; i++) {
    fprintf(output_file, "atomX %g %g %g\n", r[i].x, r[i].y, r[i].z);
  }
}


void write_vi(FILE *output_file) {
  fprintf(output_file, "%g %g %g\n", vi.x, vi.y, vi.z);
}


void write_durata_totale(void) {
  fprintf(totdurationfile,"%g\n%i\n",t*dt + last_durata_totale, nrun);
}
