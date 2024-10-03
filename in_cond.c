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


void fcc(void) {
  if (npart != 4*nlayers*npartx*nparty) {
    printf("** PROBLEMA **\nnpart non congruente!\n");
    fprintf(logfile,"** PROBLEMA **\nnpart non congruente!\n");
  }
  vec r, v0; // giochetto velocità opposte a coppie (mom angolare e mom lineare circa nulli)
  float vxm = 0., vym = 0., vzm = 0.;
  char AtomName[] = "atomX";
  float semi_lattice_len_x =  npartx*a_lattice/2.; // (0.5 + max(npartx-1,nparty-1,nlayers-1) * a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  float semi_lattice_len_y =  nparty*a_lattice/2.;
  float semi_lattice_len_z = nlayers*a_lattice/2.;
  
  cond_in = fopen("./data/in_cond.xyz","w");
  cond_in_vel = fopen("./data/in_cond_vel.dat","w");

  fprintf(cond_in,"%i\n\n",npart);

  /**** To create FCC 100 lattice*********/
  for(int i = 0; i < npartx; i++){            // Number of Atoms in the X direction

    for(int j = 0; j < nparty; j++){        // Number of Atoms in the Y direction

      for(int k = 0; k < nlayers; k++) {       // Number of layers in the Z direction

        r.x = i * a_lattice - semi_lattice_len_x;   r.y = j * a_lattice - semi_lattice_len_y;   r.z = k * a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);

        r.x = i * a_lattice - semi_lattice_len_x;   r.y = 0.5 * a_lattice +  j * a_lattice - semi_lattice_len_y;    r.z = 0.5 * a_lattice + k * a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);

        r.x = 0.5 * a_lattice + i * a_lattice - semi_lattice_len_x;  r.y = j * a_lattice - semi_lattice_len_y;   r.z = 0.5 * a_lattice + k *a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);

        r.x = 0.5 * a_lattice + i*a_lattice - semi_lattice_len_x;  r.y = 0.5 * a_lattice + j * a_lattice - semi_lattice_len_y;   r.z = k * a_lattice - semi_lattice_len_z;
        fprintf(cond_in,"%s %lf %lf %lf\n", AtomName, r.x, r.y, r.z);
      }

     }

  }

  srand(time(NULL));
  for (int i = 0; i < npart / 2; i++) { // npart / 2 perchè per annullare il momento e il momento angolare totale metto le velocità opposte a coppie
    v0.x = GaussianNoise(mu, var);
    v0.y = GaussianNoise(mu, var);
    v0.z = GaussianNoise(mu, var);
    fprintf(cond_in_vel,"%g %g %g\n",v0.x,v0.y,v0.z); // giochetto velocità opposte a coppie
    fprintf(cond_in_vel,"%g %g %g\n",-v0.x,-v0.y,-v0.z);
  }
  fclose(cond_in);
  fclose(cond_in_vel);
}


void load_r(vec *r) {
  inputr = fopen("./data/in_cond.xyz","r");
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
  fprintf(stat,"%g %g %g %g %g %g %g %g\n", t*dt + last_durata_totale, sumv.x, sumv.y, sumv.z, kenergy, penergy, kenergy + penergy, reduced_temperature);
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
  file_durata_totale = fopen("./data/durata_totale.dat","w");
  fprintf(file_durata_totale,"%g\n%i\n",t*dt + last_durata_totale, nrun);
  fclose(file_durata_totale);
}
