#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

// Il concetto è che se rmax < BOXL2, dove rmax è il troncamento del potenziale,
// allora ogni particella A interagirà solamente con una delle infinite
// immagini della stessa particella B (dimostrato con la disugluaglianza triangolare).
// Per rendere tutto consistente il reticolo creato deve riempire il box e non andare oltre.


float potenergy(float r, float eps, float sigma) {
  return 4.0*eps*(pow(sigma/r,6.) - pow(sigma/r,3.));
  // return 0.0; // per test
}


float simforce(float r, float eps, float sigma) {
  return 24.0*eps*(2.*pow(sigma/r,7.) - pow(sigma/r,4.)); // c'è l'r^2 a dividere già dentro
  // return 0.0; // per test
}

// VERLET
void verlet_periodic(vec *r, vec *ro, vec *a, params p) {
  vec rni;
  for (int i = 0; i < p.npart; i++) { // integrazione quando non devo scrivere nel file
    rni.x = 2.0 * r[i].x - ro[i].x + p.dtsquare * a[i].x / p.m;
    rni.y = 2.0 * r[i].y - ro[i].y + p.dtsquare * a[i].y / p.m;
    rni.z = 2.0 * r[i].z - ro[i].z + p.dtsquare * a[i].z / p.m;
          // printf("VERLET%g %g %g\n",rni.x,rni.y,rni.z);
    ro[i].x = r[i].x;
    ro[i].y = r[i].y;
    ro[i].z = r[i].z;
    r[i].x = rni.x - p.BOXL * round( rni.x / p.BOXL );
    r[i].y = rni.y - p.BOXL * round( rni.y / p.BOXL );
    r[i].z = rni.z - p.BOXL * round( rni.z / p.BOXL );
    a[i].x = 0.;
    a[i].y = 0.;
    a[i].z = 0.;
  }
}


void compute_kenergy_momentum(float t, vec *r, vec *ro, vec *a, float penergy, params p) {
  vec vi, sumv;
  float kenergy;
  sumv.x = 0.;
  sumv.y = 0.;
  sumv.z = 0.;
  kenergy = 0.;
  for (int i = 0; i < p.npart; i++) {
    vi.x = (r[i].x-ro[i].x-p.BOXL * round((r[i].x-ro[i].x)/p.BOXL))/p.dt; // faccio anche qui round per togliere velocità da un lato del box all'altro.
    vi.y = (r[i].y-ro[i].y-p.BOXL * round((r[i].y-ro[i].y)/p.BOXL))/p.dt;
    vi.z = (r[i].z-ro[i].z-p.BOXL * round((r[i].z-ro[i].z)/p.BOXL))/p.dt;
    sumv.x += vi.x;
    sumv.y += vi.y;
    sumv.z += vi.z;
    kenergy += vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
  }
  write_stat(t, sumv, kenergy, penergy, p);
}


void compute_kenergy_momentum2(float t, vec *r, vec *v, float penergy, params p) {
  vec sumv;
  float kenergy;
  sumv.x = 0.;
  sumv.y = 0.;
  sumv.z = 0.;
  kenergy = 0.;
  for (int i = 0; i < p.npart; i++) {
    sumv.x += v[i].x;
    sumv.y += v[i].y;
    sumv.z += v[i].z;
    kenergy += v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z;
  }
  write_stat(t, sumv, kenergy, penergy, p);
}


void compute_forces(vec *r, vec *a, params p) {
  vec rij;
  float rij2, simforceij;
  for (int i = 0; i < p.npart; i++) {
    a[i].x = 0.;
    a[i].y = 0.;
    a[i].z = 0.;
    }
  for (int i = 0; i < p.npart; i++) {
    for (int j = i + 1; j < p.npart; j++) {
      rij.x = (r[i].x - r[j].x);
      rij.y = (r[i].y - r[j].y);
      rij.z = (r[i].z - r[j].z);
      // printf("FORCES DISTANCES %g %g %g\n",rij.x,rij.y,rij.z);
      rij.x = rij.x - p.BOXL*round(rij.x/p.BOXL);
      rij.y = rij.y - p.BOXL*round(rij.y/p.BOXL);
      rij.z = rij.z - p.BOXL*round(rij.z/p.BOXL);
      rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
      // printf("FORCES DISTANCES RESCALED %g %g %g\n",rij.x,rij.y,rij.z);
      // printf("FORCES DISTANCES RESCALED MOD %g\n",rij2);
      if (rij2 < p.r_max_squared) {
        simforceij = simforce(rij2, p.eps, p.sigma);
        // printf("SIMFORCE %g\n",simforceij);
        a[i].x = a[i].x + simforceij * rij.x;
        a[i].y = a[i].y + simforceij * rij.y;
        a[i].z = a[i].z + simforceij * rij.z;
        a[j].x = a[j].x - simforceij * rij.x;
        a[j].y = a[j].y - simforceij * rij.y;
        a[j].z = a[j].z - simforceij * rij.z;
      }
    }
    // printf("FORCES %g %g %g\n",a[i].x,a[i].y,a[i].z);
  }
}


float compute_forces_stat(vec *r, vec *a, params p) {
  vec rij;
  float rij2, simforceij;
  float penergy = 0.;
  for (int i = 0; i < p.npart; i++) {
    a[i].x = 0.;
    a[i].y = 0.;
    a[i].z = 0.;
    }
  for (int i = 0; i < p.npart; i++) { // calcolo forze
    for (int j = i + 1; j < p.npart; j++) {
      rij.x = (r[i].x - r[j].x);
      rij.y = (r[i].y - r[j].y);
      rij.z = (r[i].z - r[j].z);
      rij.x = rij.x - p.BOXL*round(rij.x/p.BOXL);
      rij.y = rij.y - p.BOXL*round(rij.y/p.BOXL);
      rij.z = rij.z - p.BOXL*round(rij.z/p.BOXL);
      rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
      if (rij2 < p.r_max_squared) {
        simforceij = simforce(rij2, p.eps, p.sigma);
        penergy = penergy + (potenergy(rij2, p.eps, p.sigma) - p.shift);
        a[i].x = a[i].x + simforceij * rij.x;
        a[i].y = a[i].y + simforceij * rij.y;
        a[i].z = a[i].z + simforceij * rij.z;
        a[j].x = a[j].x - simforceij * rij.x;
        a[j].y = a[j].y - simforceij * rij.y;
        a[j].z = a[j].z - simforceij * rij.z;
      }
    }
  }
  return penergy;
}


void eulero(vec *r, vec *ro, vec *a, params p) {
  vec vi, rni;
  FILE *inputv = fopen(p.restartvelfilepath,"r");
  for (int i = 0; i < p.npart; i++) {
    fscanf(inputv, "%g %g %g\n", &(vi.x), &(vi.y), &(vi.z));
    rni.x = r[i].x + vi.x * p.dt + 0.5 * p.dtsquare * a[i].x / p.m;
    rni.y = r[i].y + vi.x * p.dt + 0.5 * p.dtsquare * a[i].y / p.m;
    rni.z = r[i].z + vi.x * p.dt + 0.5 * p.dtsquare * a[i].z / p.m;
    // printf("EULERO %g %g %g\n",rni.x,rni.y,rni.z);
    ro[i].x = r[i].x;
    ro[i].y = r[i].y;
    ro[i].z = r[i].z;
    r[i].x = rni.x - p.BOXL * round(rni.x/p.BOXL);
    r[i].y = rni.y - p.BOXL * round(rni.y/p.BOXL);
    r[i].z = rni.z - p.BOXL * round(rni.y/p.BOXL);
  }
  fclose(inputv);
}
