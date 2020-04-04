#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"


int main() {
	// APRO FILES E LI INIZIALIZZO
	param = fopen("./param.in","r");
	inputr = fopen("./data/verlet_periodic.xyz","r");
	output = fopen("./data/corr_func.dat","w");

	// CARICO PARAMETRI DA PARAM.IN
	fscanf(param,"npartx=%i\nnparty=%i\nnlayers=%i\nnpart=%i\nwrite_jump=%i\ntimesteps=%i\ndt=%g\neps=%g\nsigma=%g\nmu=%g\nvar=%g\nm=%g\na_lattice=%g\npot_trunc_perc=%g\nnew_in_cond=%i",&npartx,&nparty,&nlayers,&npart,&write_jump,&timesteps,&dt,&eps,&sigma,&mu,&var,&m,&a_lattice,&pot_trunc_perc,&newc);
	fclose(param);

	// INIZIALIZZO VARIABILI DA PARAMETRI
	dtdouble = 2.*dt;
	dtsquare = pow(dt, 2.);
	r_max = sigma * pow((1+sqrt( 1-16*pot_trunc_perc ))/( 2*pot_trunc_perc ), 1./6.);	/*impostando la percentuale di troncamento del potenziale, posso calcolare rmax sapendo che è il valore in cui il potenziale raggiunge quella percentuale impostata */
	r_max_squared = pow(r_max, 2.);
	shift = potenergy(r_max, eps, sigma); // per lo shift del potenziale
	BOXL = ( 0.5 + nlayers ) * a_lattice; // (0.5 + max(npartx,nparty,nlayers) * a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
	reduced_density = npart * sigma / pow(BOXL, 3.);









	// CARICO DURATA TOTALE
	file_durata_totale = fopen("./data/durata_totale.dat","r");
	fscanf(file_durata_totale,"%g\n%i\n",&last_durata_totale,&nrun);
	fclose(file_durata_totale);
	
	// INIZIALIZZO VARIABILI AUTODIFFUSION
  	int ntime = timesteps / write_jump;
	float dx = sigma / 30.;
	int n_bin = (int)(r_max / dx);
	float xv[n_bin], fv[n_bin];
	for (int a = 0; a < n_bin; a++) xv[a] = dx*a + 0.5*dx;
	for (int a = 0; a < n_bin; a++) fv[a] = 0;
	vec r[npart], rij;
	float rij2;

	// INIZIALIZZO FILE OUTPUT
	fprintf(output, "#    x(1)      f(x)(2)\n");

    // CARICO r
	for (int k = 0; k < ntime; k++) {
		load_r_2(inputr, r);
		
		// CICLO SU OGNI PARTICELLA PRIMARIA j
		for (int j = 0; j < npart; j++) {
			
			// CICLO SU OGNI PARTICELLA SECONDARIA i
			for (int i = 0; i < npart; i++) {
				rij.x = (r[j].x - r[i].x) - BOXL*round((r[j].x-r[i].x)/BOXL);
				rij.y = (r[j].y - r[i].y) - BOXL*round((r[j].y-r[i].y)/BOXL);
				rij.z = (r[j].z - r[i].z) - BOXL*round((r[j].z-r[i].z)/BOXL);
				rij2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
				if (rij2 < r_max_squared && rij2 != 0.) {
					rij2 = pow(rij2, 0.5);
					fv[(int)floor(rij2/dx)]++;
				}
			}
		}
	}
	for (int a = 0; a < n_bin; a++) {
		fv[a] /= (npart*ntime);
		fprintf(output, "%g %g %g\n",xv[a],fv[a],pow(xv[a],2.)*4.*3.1415*dx);
	}
	fclose(output);
}

