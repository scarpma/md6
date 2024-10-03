#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"


int main() {
  // APRO FILES E LI INIZIALIZZO
  param = fopen("./param.in","r");
  inputr = fopen("./data/verlet_periodic.xyz","r");
  output = fopen("./data/autodiffusion.dat","w");

  // CARICO PARAMETRI DA PARAM.IN
  fscanf(param,"npartx=%i\nnparty=%i\nnlayers=%i\nnpart=%i\nwrite_jump=%i\ntimesteps=%i\ndt=%g\neps=%g\nsigma=%g\nmu=%g\nvar=%g\nm=%g\na_lattice=%g\npot_trunc_perc=%g\nnew_in_cond=%i",&npartx,&nparty,&nlayers,&npart,&write_jump,&timesteps,&dt,&eps,&sigma,&mu,&var,&m,&a_lattice,&pot_trunc_perc,&newc);
  fclose(param);

  // INIZIALIZZO VARIABILI DA PARAMETRI
  dtdouble = 2.*dt;
  dtsquare = pow(dt, 2.);
  r_max = sigma * pow((1+sqrt( 1-16*pot_trunc_perc ))/( 2*pot_trunc_perc ), 1./6.);  /*impostando la percentuale di troncamento del potenziale, posso calcolare rmax sapendo che è il valore in cui il potenziale raggiunge quella percentuale impostata */
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
  vec r[npart], rin[npart], disp2[npart], var;
  float time;

  // INIZIALIZZO FILE OUTPUT
  fprintf(output, "#    t(1)      varx(2)      vary(3)      varz(3)\n");
  
  // CICLO SUI TEMPI

  for (int j = 0; j < ntime; j++) { // ciclo su ogni istante di tempo (scritto nel file)

    // CARICO r e rin
    load_r_2(inputr, r);
    // for (int i = 0; i < npart; i++) printf("%g %g %g\n",r[i].x,r[i].y,r[i].z);
    if (j == 0) {  
      for (int i = 0; i < npart; i++) {
        rin[i].x = r[i].x;
        rin[i].y = r[i].y;
        rin[i].z = r[i].z;
      }
    }


    // CICLO SULLE PARTICELLE AL TEMPO j (NO ! Non va bene così. Perchè comunque rx -rxin non sarà mai maggiore di ????)
    for (int i = 0; i < npart; i++) {
      disp2[i].x = pow(( r[i].x - rin[i].x ) - BOXL*round(( r[i].x - rin[i].x )/BOXL), 2.);
      disp2[i].y = pow(( r[i].y - rin[i].y ) - BOXL*round(( r[i].y - rin[i].y )/BOXL), 2.);
      disp2[i].z = pow(( r[i].z - rin[i].z ) - BOXL*round(( r[i].z - rin[i].z )/BOXL), 2.);
      //printf("%g %g %g\n",disp2[i].x,disp2[i].y,disp2[i].z);          
    }
    var.x = 0.; var.y = 0.; var.z = 0.;
    for (int i = 0; i < npart; i++) {
      var.x += disp2[i].x;
      var.y += disp2[i].y;
      var.z += disp2[i].z;
    }
    time = j * dt * write_jump + last_durata_totale - timesteps*dt;
    var.x = var.x/npart;
    var.y = var.y/npart;
    var.z = var.z/npart;
    fprintf(output, "%g %g %g %g\n",time,var.x,var.y,var.z);

  }
}





