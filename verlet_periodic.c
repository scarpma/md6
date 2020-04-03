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
void verlet_periodic(vec *r, vec *ro, vec *a) {
	for (int i = 0; i < npart; i++) { // integrazione quando non devo scrivere nel file
			rni.x = 2.0 * r[i].x - ro[i].x + dtsquare * a[i].x / m;
			rni.y = 2.0 * r[i].y - ro[i].y + dtsquare * a[i].y / m;
			rni.z = 2.0 * r[i].z - ro[i].z + dtsquare * a[i].z / m;
			ro[i].x = r[i].x;
			ro[i].y = r[i].y;
			ro[i].z = r[i].z;
			r[i].x = rni.x - BOXL * round( rni.x / BOXL );
			r[i].y = rni.y - BOXL * round( rni.y / BOXL );
			r[i].z = rni.z - BOXL * round( rni.z / BOXL );
			a[i].x = 0.;
			a[i].y = 0.;
			a[i].z = 0.;
		}
}

// VERLET CON CALCOLO MOMENTO ED ENERGIA TOTALE
void verlet_periodic_write(vec *r, vec *ro, vec *a) {
	for (int i = 0; i < npart; i++) {
			rni.x = 2.0 * r[i].x - ro[i].x + dtsquare * a[i].x / m;
			rni.y = 2.0 * r[i].y - ro[i].y + dtsquare * a[i].y / m;
			rni.z = 2.0 * r[i].z - ro[i].z + dtsquare * a[i].z / m;
			vi.x = (rni.x-ro[i].x-BOXL * round((rni.x-ro[i].x)/BOXL))/dtdouble; // faccio anche qui round per togliere velocità da un lato del box all'altro.
			vi.y = (rni.y-ro[i].y-BOXL * round((rni.y-ro[i].y)/BOXL))/dtdouble;
			vi.z = (rni.z-ro[i].z-BOXL * round((rni.z-ro[i].z)/BOXL))/dtdouble;
			ro[i].x = r[i].x;
			ro[i].y = r[i].y;
			ro[i].z = r[i].z;
			r[i].x = rni.x - BOXL * round( rni.x / BOXL );
			r[i].y = rni.y - BOXL * round( rni.y / BOXL );
			r[i].z = rni.z - BOXL * round( rni.z / BOXL );
			a[i].x = 0.;
			a[i].y = 0.;
			a[i].z = 0.;
			sumv.x += vi.x;
			sumv.y += vi.y;
			sumv.z += vi.z;
			kenergy += vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
		}
	write_r(output, r);
	write_stat();
	sumv.x = 0.0;
	sumv.y = 0.0;
	sumv.z = 0.0;
	penergy = 0.0;
	kenergy = 0.0;
}


void verlet_periodic_last(vec *r, vec *ro, vec *a) {
	for (int i = 0; i < npart; i++) {
		rni.x = 2.0 * r[i].x - ro[i].x + dtsquare * a[i].x / m;
		rni.y = 2.0 * r[i].y - ro[i].y + dtsquare * a[i].y / m;
		rni.z = 2.0 * r[i].z - ro[i].z + dtsquare * a[i].z / m;
		vi.x = (rni.x-ro[i].x-BOXL * round((rni.x-ro[i].x)/BOXL))/dtdouble; // faccio anche qui round per togliere velocità da un lato del box all'altro.
		vi.y = (rni.y-ro[i].y-BOXL * round((rni.y-ro[i].y)/BOXL))/dtdouble;
		vi.z = (rni.z-ro[i].z-BOXL * round((rni.z-ro[i].z)/BOXL))/dtdouble;
		ro[i].x = r[i].x;
		ro[i].y = r[i].y;
		ro[i].z = r[i].z;
		r[i].x = rni.x - BOXL * round( rni.x / BOXL );
		r[i].y = rni.y - BOXL * round( rni.y / BOXL );
		r[i].z = rni.z - BOXL * round( rni.z / BOXL );
		a[i].x = 0.;
		a[i].y = 0.;
		a[i].z = 0.;
		sumv.x += vi.x;
		sumv.y += vi.y;
		sumv.z += vi.z;
		kenergy += vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
		write_vi(output_in_cond_vel);
	}
	write_r(output, r);
	write_r(output_in_cond, r);
	write_stat();
}


void compute_forces(vec *r, vec *a) {
	for (int i = 0; i < npart; i++) {
		for (int j = i + 1; j < npart; j++) {
			rij.x = (r[i].x - r[j].x);
			rij.y = (r[i].y - r[j].y);
			rij.z = (r[i].z - r[j].z);
			rij.x = rij.x - BOXL*round(rij.x/BOXL);
			rij.y = rij.y - BOXL*round(rij.y/BOXL);
			rij.z = rij.z - BOXL*round(rij.z/BOXL);
			rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
			if (rij2 < r_max_squared) {
				simforceij = simforce(rij2, eps, sigma);
				a[i].x = a[i].x + simforceij * rij.x;
				a[i].y = a[i].y + simforceij * rij.y;
				a[i].z = a[i].z + simforceij * rij.z;
				a[j].x = a[j].x - simforceij * rij.x;
				a[j].y = a[j].y - simforceij * rij.y;
				a[j].z = a[j].z - simforceij * rij.z;
			}
		}
	}
}


void compute_forces_stat(vec *r, vec *a) {
	for (int i = 0; i < npart; i++) { // calcolo forze
		for (int j = i + 1; j < npart; j++) {
			rij.x = (r[i].x - r[j].x);
			rij.y = (r[i].y - r[j].y);
			rij.z = (r[i].z - r[j].z);
			rij.x = rij.x - BOXL*round(rij.x/BOXL);
			rij.y = rij.y - BOXL*round(rij.y/BOXL);
			rij.z = rij.z - BOXL*round(rij.z/BOXL);
			rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
			if (rij2 < r_max_squared) {
				simforceij = simforce(rij2, eps, sigma);
				penergy = penergy + (potenergy(rij2, eps, sigma) - shift);
				a[i].x = a[i].x + simforceij * rij.x;
				a[i].y = a[i].y + simforceij * rij.y;
				a[i].z = a[i].z + simforceij * rij.z;
				a[j].x = a[j].x - simforceij * rij.x;
				a[j].y = a[j].y - simforceij * rij.y;
				a[j].z = a[j].z - simforceij * rij.z;
			}
		}
	}
}


void eulero(vec *r, vec *ro, vec *a) {
	vec vi, rni;
	inputv = fopen("./data/in_cond_vel.dat","r");
	for (int i = 0; i < npart; i++) {
		fscanf(inputv, "%g %g %g\n", &(vi.x), &(vi.y), &(vi.z));
		rni.x = r[i].x + vi.x * dt + 0.5 * dtsquare * a[i].x / m;
		rni.y = r[i].y + vi.x * dt + 0.5 * dtsquare * a[i].y / m;
		rni.z = r[i].z + vi.x * dt + 0.5 * dtsquare * a[i].z / m;
		ro[i].x = r[i].x;
		ro[i].y = r[i].y;
		ro[i].z = r[i].z;
		r[i].x = rni.x - BOXL * round(rni.x/BOXL);
		r[i].y = rni.y - BOXL * round(rni.y/BOXL);
		r[i].z = rni.z - BOXL * round(rni.y/BOXL);
		a[i].x = 0.;
		a[i].y = 0.;
		a[i].z = 0.;
	}
	fclose(inputv);
}
