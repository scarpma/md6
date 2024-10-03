#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

float potenergy(float r, float eps, float sigma) {
  return 4.0*eps*(pow(sigma/r,6.) - pow(sigma/r,3.));
  // return 0.0; // per test
}

typedef struct {
  float x,y,z;
} vec;

float GaussianNoise(float mu, float var){
  double u1, u2;
  u1 = rand() * (1.0 / RAND_MAX);
  u2 = rand() * (1.0 / RAND_MAX);
  double z0, z1;
  z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(2 * M_PI * u2);
  return z0 * var + mu;
}

void fcc(int nlayers, int npartx, int nparty, float a_lattice, int npart, float mu, float var) {
  if (npart != 4*nlayers*npartx*nparty) {
    printf("** PROBLEMA **\nnpart non congruente!\n");
  }
  vec r, v0; // giochetto velocità opposte a coppie (mom angolare e mom lineare circa nulli)
  float vxm = 0., vym = 0., vzm = 0.;
  char AtomName[] = "atomX";
  float semi_lattice_len_x =  npartx*a_lattice/2.; // (0.5 + max(npartx-1,nparty-1,nlayers-1) * a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  float semi_lattice_len_y =  nparty*a_lattice/2.;
  float semi_lattice_len_z = nlayers*a_lattice/2.;
  
  FILE *cond_in;
  cond_in = fopen("./data/in_cond.xyz","w");
  FILE *cond_in_vel;
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

int main() {

  int npartx, nparty, nlayers, npart, write_jump, newc;
  int timesteps;
  float dt, eps, sigma, mu, var, m; 
  float a_lattice,pot_trunc_perc;
  float r_max, BOXL, reduced_density;
  float dtdouble, dtsquare, r_max_squared, shift;
  

  // APRO FILE DEI PARAMETRI
  FILE *param;
  param = fopen("./param.in","r");

  // CARICO PARAMETRI DA PARAM.IN
  fscanf(param, 
    "npartx=%i\n"
    "nparty=%i\n"
    "nlayers=%i\n"
    "npart=%i\n"
    "write_jump=%i\n"
    "timesteps=%i\n"
    "dt=%g\n"
    "eps=%g\n"
    "sigma=%g\n"
    "mu=%g\n"
    "var=%g\n"
    "m=%g\n"
    "a_lattice=%g\n"
    "pot_trunc_perc=%g\n"
    "new_in_cond=%i",
    &npartx, &nparty, &nlayers, &npart, &write_jump, 
    &timesteps, &dt, &eps, &sigma, &mu, &var, &m, 
    &a_lattice, &pot_trunc_perc, &newc);
  fclose(param);

  // INIZIALIZZO VARIABILI DA PARAMETRI
  dtdouble = 2.*dt;
  dtsquare = pow(dt, 2.);
  r_max = sigma * pow((1+sqrt( 1-16*pot_trunc_perc ))/( 2*pot_trunc_perc ), 1./6.);  /*impostando la percentuale di troncamento del potenziale, posso calcolare rmax sapendo che è il valore in cui il potenziale raggiunge quella percentuale impostata */
  r_max_squared = pow(r_max, 2.);
  shift = potenergy(r_max, eps, sigma); // per lo shift del potenziale
  BOXL = ( 0.5 + nlayers ) * a_lattice; // (0.5 + max(npartx,nparty,nlayers) * a_lattice è la lunghezza  del reticolo nella direzione che la rende massima)
  reduced_density = npart * sigma / pow(BOXL, 3.);

  FILE *logfile;
  logfile = fopen("./salvati/run.log", "w");
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
  fclose(logfile);


  // FILE *output;
  // output = fopen("./data/verlet_periodic.xyz","w");
  vec *r, *ro, *a;
  r  = malloc(npart*sizeof(vec));
  ro = malloc(npart*sizeof(vec));
  a  = malloc(npart*sizeof(vec));

  fcc(nlayers, npartx, nparty, a_lattice, npart, mu, var);

  return 0;
}
