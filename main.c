#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"


int main() {

  // APRO FILE DEI PARAMETRI
  param = fopen("./param.in","r");
  if (!param) {
    perror("fopen param");
    exit(1);
  }

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

  // VEDO SE CONTINUARE O INIZIARE UNA NUOVA SIMULAZIONE
  if (newc == 0) {
    logfile = fopen("./salvati/run.log", "a");
    if (!logfile) {
      perror("fopen logfile");
      exit(1);
    }
    file_durata_totale = fopen("./data/durata_totale.dat","r");
    fscanf(file_durata_totale,"%g\n%i\n",&last_durata_totale,&nrun);
    nrun = nrun + 1;
    fclose(file_durata_totale);
    printf("Proseguo ultima simulazione. nrun = %d, t_in = %.2f\n",nrun,last_durata_totale);
    fprintf(logfile,"Proseguo ultima simulazione. nrun = %d, t_in = %.2f\n",nrun,last_durata_totale);
  } else if (newc == 1) {
    int status = system("./del-data"); // esegue uno script di sistema
    status = system("mkdir ./salvati"); // esegue uno script di sistema
    status = system("mkdir ./data"); // esegue uno script di sistema
    printf("r_max=%g    BOXL=%g    red. dens=%g\n",r_max,BOXL,reduced_density);
    nrun = 0;
    last_durata_totale = 0.;
    logfile = fopen("./salvati/run.log", "w");
    if (!logfile) {
      perror("fopen logfile");
      exit(1);
    }
    printf("Inizio una nuova simulazione: inizializzo reticolo FCC e velocità random.\n");
    fprintf(logfile,"Inizio una nuova simulazione:\n\nPARAM:\n");
    fprintf(logfile,"npartx=%i\nnparty=%i\n",npartx,nparty);
    fprintf(logfile,"nlayers=%i\nnpart=%i\n",nlayers,npart);
    fprintf(logfile,"write_jump=%i\ntimesteps=%i\n",write_jump,timesteps);
    fprintf(logfile,"dt=%g\neps=%g\nsigma=%g\n",dt,eps,sigma);
    fprintf(logfile,"mu=%g\nvar=%g\nm=%g\na_lattice=%g\n",mu,var,m,a_lattice);
    fprintf(logfile,"pot_trunc_perc=%g\nnew_in_cond=%i\n\n\n",pot_trunc_perc,newc);
    fprintf(logfile,"inizializzo reticolo FCC e velocità random\n\n");
    fprintf(logfile, "r_max=%g    BOXL=%g    red. dens=%g\n",r_max,BOXL,reduced_density);
    fcc();
  }
    
    
  // APRO FILE OUTPUT E LI INIZIALIZZO
  stat = fopen("./data/stat.dat","w");
  if (!stat) {
    perror("fopen stat");
    exit(1);
  }
  fprintf(stat,"# t(0)   sumvx(1)   sumvy(2)   sumvz(3)    kenergy(4)   penergy(5)   energy(6)   red_temp(7)\n");
  output = fopen("./data/verlet_periodic.xyz","w");
    
  // DICHIARO VARIABILI DINAMICHE
  vec r[npart], ro[npart], a[npart];
  
  // CARICO CONDIZIONI INIZIALI E GENERO LE ro
  load_r(r);
  write_r(output, r);
  printf("Inizio ad integrare:\n\n");
  fprintf(logfile, "Inizio ad integrare: attendi!\n\n");
    for (int i = 0; i < npart; i++) {
        a[i].x = 0.0;
        a[i].y = 0.0;
        a[i].z = 0.0;
    }
  compute_forces(r, a);
  eulero(r, ro, a);
  write_r(output, r);


  // CONTROLLO CHE TUTTE LE PARTICELLE SIANO NEL BOX
  int out_count = 0;
  for (int i = 0; i < npart; i++) {
    if (r[i].x < -BOXL/2. || r[i].x > BOXL/2.) out_count += 1;
    if (r[i].y < -BOXL/2. || r[i].y > BOXL/2.) out_count += 1;
    if (r[i].z < -BOXL/2. || r[i].z > BOXL/2.) out_count += 1;
  }
  if (out_count > 0) printf("Le condizioni iniziali non sono all'interno del box periodico.\n");


  // INTEGRO CON VERLET
  t = 2;
  for (int tt = 2; tt < (timesteps-1) / write_jump; tt++) {
    for (int k = 0; k < write_jump; k++) {
      compute_forces(r, a);
      verlet_periodic(r, ro, a);
      t++;
    }
    compute_forces_stat(r, a);
    verlet_periodic_write(r, ro, a);
    t++;
  }
  
  // ULTIMA INTEGRAZIONE
  output_in_cond = fopen("./data/in_cond.xyz","w");
    output_in_cond_vel = fopen("./data/in_cond_vel.dat","w");
  compute_forces_stat(r, a);
    verlet_periodic_last(r, ro, a); t++;
  write_durata_totale();
  
  // CHIUDO FILES
  fclose(output);
  fclose(output_in_cond);
  fclose(output_in_cond_vel);
  fclose(stat);
  printf("Finito!\n\n");
  fprintf(logfile, "Finito!\n\n\n\n\n\n\n\n");
  fclose(logfile);
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
