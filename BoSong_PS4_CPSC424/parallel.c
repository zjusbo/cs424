#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <mpi.h>
/*
  This is the serial program for CPSC424/524 Assignment #4.
  
  Author: Andrew Sherman, Yale University
  Modified by: Bo Song, Yale University
  Date: 2/23/2013

*/
// Globals so that we don't need to pass them to functions
double dt, dt2;
int N, n, K;
double *mass, *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz, *buf;
double *bmass, *bx, *by, *bz; // boundary bodies

void swap(double ** a, double ** b){
  double * tmp = *a;
  *a = *b;
  *b = tmp;
}

// copyed from serial.c
//Function to compute the forces between a pair of bodies
void force(int body1, int body2, double *deltaf) {
  double gmmr3, r, r2, dx, dy, dz;
  double G=1.0;
  
  dx = x[body2] - x[body1];
  dy = y[body2] - y[body1];
  dz = z[body2] - z[body1];
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  if (r<=5.) {
    gmmr3 = G*mass[body1]*mass[body2]/(r2*r);
    deltaf[0] = gmmr3 * dx;
    deltaf[1] = gmmr3 * dy;
    deltaf[2] = gmmr3 * dz;
  }
  else{
    deltaf[0] = 0.;
    deltaf[1] = 0.;
    deltaf[2] = 0.;
  }
} 

// modified from force()
//Function to compute the forces between a pair of body and boundary body
void bforce(int body1, int body2, double *deltaf) {
  double gmmr3, r, r2, dx, dy, dz;
  double G=1.0;
  
  dx = bx[body2] - x[body1];
  dy = by[body2] - y[body1];
  dz = bz[body2] - z[body1];
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  if (r<=5.) {
    gmmr3 = G*mass[body1]*bmass[body2]/(r2*r);
    deltaf[0] = gmmr3 * dx;
    deltaf[1] = gmmr3 * dy;
    deltaf[2] = gmmr3 * dz;
  }
  else{
    deltaf[0] = 0.;
    deltaf[1] = 0.;
    deltaf[2] = 0.;
  }
}
// Copyted from serial.c
// Function to print center of mass and average velocity
void output(int ts) {
  int thisbody;
  double cmassx, cmassy, cmassz, tmass, tvelx, tvely, tvelz;

  if (ts==0) printf("\n\nInitial Conditions (time = 0.0):\n");
  else printf("\n\nConditions after timestep %d (time = %f):\n",ts,ts*dt);

  cmassx = 0.;
  cmassy = 0.;
  cmassz = 0.;
  tmass = 0.;
  tvelx = 0.;
  tvely = 0.;
  tvelz = 0.;
  for (thisbody=0; thisbody<n; thisbody++) {
    cmassx += mass[thisbody]*x[thisbody];
    cmassy += mass[thisbody]*y[thisbody];
    cmassz += mass[thisbody]*z[thisbody];
    tmass  += mass[thisbody];
    tvelx += vx[thisbody];
    tvely += vy[thisbody];
    tvelz += vz[thisbody];
  }
  printf("\n     Center of Mass:   (%e, %e, %e)\n", cmassx/tmass, cmassy/tmass, cmassz/tmass);
  printf(  "     Average Velocity: (%e, %e, %e)\n", tvelx/n, tvely/n, tvelz/n);
}

void intermediate_output(int ts, int rank){
  int thisbody;
  double cmassx, cmassy, cmassz, tmass, tvelx, tvely, tvelz; 
  int i;
  if(rank == 0){
    if (ts==0) printf("\n\nInitial Conditions (time = 0.0):\n");
    else printf("\n\nConditions after timestep %d (time = %f):\n",ts,ts*dt);
  }
  cmassx = 0.;
  cmassy = 0.;
  cmassz = 0.;
  tmass = 0.;
  tvelx = 0.;
  tvely = 0.;
  tvelz = 0.;
  for (thisbody=0; thisbody<n; thisbody++) {
    cmassx += mass[thisbody]*x[thisbody];
    cmassy += mass[thisbody]*y[thisbody];
    cmassz += mass[thisbody]*z[thisbody];
    tmass  += mass[thisbody];
    tvelx += vx[thisbody];
    tvely += vy[thisbody];
    tvelz += vz[thisbody];
  }
  int tmp[8];
  MPI_Gather(&n, 1, MPI_INT, tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(rank == 0){
    MPI_Reduce(MPI_IN_PLACE, &cmassx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &cmassy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &cmassz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &tmass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &tvelx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &tvely, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &tvelz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }else{
    MPI_Reduce(&cmassx, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cmassy, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cmassz, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tmass, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tvelx, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tvely, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tvelz, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  }
  if(rank == 0){
    printf("\n     Center of Mass:   (%e, %e, %e)\n", cmassx/tmass, cmassy/tmass, cmassz/tmass);
    printf(  "     Average Velocity: (%e, %e, %e)\n", tvelx/N, tvely/N, tvelz/N);
    printf("	Load balance:");
    for(i = 0; i < 8; i++){
      printf("%d ", tmp[i]);
    }
    printf("\n");
  }
}

// partition bodies into 8 groups based on its 3D coordinates. 
// offsets is an int array containing 8 elements, each of them points to the start position of a group
// x, y, z, vx, vy, vz are global variables. 
// they can be accessed in this function
void partition(int n, int* offsets){
  // use fx, fy, fz, bx, by, bz, bmass as buffers
  int i, j;
  int sector[][3] = {{1,1,1}, {1,1,-1}, {1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
  int p = 0;
  for(i = 0; i < 8; i++){
    offsets[i] = p;
    for(j = 0; j < n; j++){
      // potential bug here
      if(x[j] * sector[i][0] > 0 && y[j] * sector[i][1] > 0 && z[j] * sector[i][2] > 0){
        bx[p] = x[j];
        by[p] = y[j];
        bz[p] = z[j];
        fx[p] = vx[j];
        fy[p] = vy[j];
        fz[p] = vz[j];
        bmass[p] = mass[j];
        p++;
      }
    }
  }
  swap(&x, &bx);
  swap(&y, &by);
  swap(&z, &bz);
  swap(&vx, &fx);
  swap(&vy, &fy);
  swap(&vz, &fz);
  swap(&mass, &bmass);
}

// extract and copy close-to-boundary points to new arrays
// bx, by, bz, bmass
void extract(int n, int* offsets, int* bn){
  int sector[][3] = {{1,1,1}, {1,1,-1}, {1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
  int i, j;
  int p = 0;
  for(i = 0; i < 8; i++){
    offsets[i] = p;
    for(j = 0; j < n; j++){
      if(x[j] * sector[i][0] > 0 && y[j] * sector[i][1] > 0 && z[j] * sector[i][2] > 0){
        // body is in current cube
        continue;
      }
      if((x[j] * x[j] + y[j] * y[j] + z[j] * z[j] <= 25.) ||
          (x[j] * x[j] + y[j] * y[j] <= 25. && z[j] * sector[i][2] > 0) ||
          (x[j] * x[j] + z[j] * z[j] <= 25. && y[j] * sector[i][1] > 0) ||
          (z[j] * z[j] + y[j] * y[j] <= 25. && x[j] * sector[i][0] > 0) ||
          (x[j] * x[j] <= 25. && y[j] * sector[i][1] > 0 && z[j] * sector[i][2] > 0) ||
          (y[j] * y[j] <= 25. && x[j] * sector[i][0] > 0 && z[j] * sector[i][2] > 0) ||
          (z[j] * z[j] <= 25. && x[j] * sector[i][0] > 0 && y[j] * sector[i][1] > 0)
        ){
        // add current body to boundary body set
        bx[p] = x[j];
        by[p] = y[j];
        bz[p] = z[j];
        bmass[p] = mass[j];
        p++;
      }
    }
  }
  *bn = p;
}

// Main function
int main(int argc, char **argv) {
  int i, ts, thisbody, otherbody;
  int bn;
  double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
  // Serial Timing Functions: See the timer.h include file in this directory
  double etime, etime0, etime1, cptime;

  int send_counts[8], recv_counts[8], send_offsets[8], recv_offsets[8];  
  int rank, size, type=99;
  MPI_Init(&argc, &argv);  
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Get no. of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?
  
  if(rank == 0){
    // Read basic inputs   
    scanf("%d\n",&N);
    scanf("%d\n",&K);
    scanf("%le\n",&dt);
  }
  
  // broadcast N, K, dt to all processes
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  dt2 = dt/2.;

  // Allocate arrays
  // times 8 to make sure it won't excced the boundary use. 
  // since one body may be boundary body for 7 processes
  mass = calloc(N * 8, sizeof(double)); // mass[i] is mass of body i
  x = calloc(N * 8, sizeof(double)); // x[i] is x position of body i
  y = calloc(N * 8, sizeof(double)); // y[i] is y position of body i
  z = calloc(N * 8, sizeof(double)); // z[i] is z position of body i
  vx = calloc(N * 8, sizeof(double)); // vx[i] is x velocity of body i
  vy = calloc(N * 8, sizeof(double)); // vy[i] is y velocity of body i
  vz = calloc(N * 8, sizeof(double)); // vz[i] is z velocity of body i
  fx = calloc(N * 8, sizeof(double)); // fx[i] is x force on body i
  fy = calloc(N * 8, sizeof(double)); // fy[i] is y force on body i
  fz = calloc(N * 8, sizeof(double)); // fz[i] is z force on body i

  bmass = calloc(N * 8, sizeof(double)); // boundary body mass
  bx = calloc(N * 8, sizeof(double)); // boundary body x
  by = calloc(N * 8, sizeof(double)); // boundary body y
  bz = calloc(N * 8, sizeof(double)); // boundary body z
  
  buf = calloc(N * 8, sizeof(double)); // boundary body z

  if(rank == 0){
    // Read initial conditions
    for (i=0; i<N; i++) scanf("%le\n",&mass[i]);
    for (i=0; i<N; i++) scanf("%le %le %le\n",&x[i],&y[i],&z[i]);
    for (i=0; i<N; i++) scanf("%le %le %le\n",&vx[i],&vy[i],&vz[i]);
    n = N;
  }else{
    n = 0;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  // Start Timer
  timing(&etime0,&cptime);

  // Timestamp Loop
  // loop phase 
  for (ts=0; ts<K; ts++) {
    if (ts%128 == 0) intermediate_output(ts, rank); // Print output if necessary
    // partition bodies into 8 groups based on its coordination
    partition(n, send_offsets);
    for(i = 0; i < size - 1; i++){
      send_counts[i] = send_offsets[i + 1] - send_offsets[i];
    }
    send_counts[size - 1] = n - send_offsets[size - 1];
    // alltoall number of bodies
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
    recv_offsets[0] = 0;
    for(i = 1; i < size; i++){
      recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
    }
    // n is not correct
    // update n
    n = recv_offsets[size - 1] + recv_counts[size - 1];
    // alltoallv bodies
    MPI_Alltoallv(mass, send_counts, send_offsets, MPI_DOUBLE, 
      buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&mass, &buf);
    MPI_Alltoallv(x, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&x, &buf);
    MPI_Alltoallv(y, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&y, &buf);
    MPI_Alltoallv(z, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&z, &buf);
    MPI_Alltoallv(vx, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&vx, &buf);
    MPI_Alltoallv(vy, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&vy, &buf);
    MPI_Alltoallv(vz, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&vz, &buf);
    // extract boundary points
    extract(n, send_offsets, &bn);
   // alltoall number of boundary points
    for(i = 0; i < size - 1; i++){
      send_counts[i] = send_offsets[i + 1] - send_offsets[i];
    }
    send_counts[size - 1] = bn - send_offsets[size - 1];
    // alltoall number of bodies
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
     
    recv_offsets[0] = 0;
    for(i = 1; i < size; i++){
      recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
    }
    bn = recv_offsets[size - 1] + recv_counts[size - 1];
    
    // alltoallv boundary points mass
    MPI_Alltoallv(bmass, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&bmass, &buf);

    // alltoallv boundary points x
    MPI_Alltoallv(bx, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&bx, &buf);
    
    // alltoallv boundary points y
    MPI_Alltoallv(by, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&by, &buf);

    // alltoallv boundary points z
    MPI_Alltoallv(bz, send_counts, send_offsets, MPI_DOUBLE, 
          buf, recv_counts, recv_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    swap(&bz, &buf);
        // Initialize forces on bodies
    for (thisbody=0; thisbody<n; thisbody++) {
      fx[thisbody] = 0.;
      fy[thisbody] = 0.;
      fz[thisbody] = 0.;
    }

    // Compute all pairwise interbody forces (Note that we take advantage of symmetry.)
    for (thisbody=0; thisbody<n; thisbody++) {
      for (otherbody=thisbody+1; otherbody<n; otherbody++) {
        force(thisbody, otherbody, deltaf); // This function computes the pairwise force of otherbody on thisbody
        fx[thisbody] += deltaf[0]; // Add x component of force to thisbody
        fy[thisbody] += deltaf[1]; // Add y component of force to thisbody
        fz[thisbody] += deltaf[2]; // Add z component of force to thisbody
        fx[otherbody] -= deltaf[0]; // Subtract x component of force from otherbody
        fy[otherbody] -= deltaf[1]; // Subtract y component of force from otherbody
        fz[otherbody] -= deltaf[2]; // Subtract z component of force from otherbody
      }
    }

    // Compute all pairwise boundary body forces (Note that we take advantage of symmetry.)
    for (thisbody=0; thisbody<n; thisbody++) {
      for (otherbody=0; otherbody<bn; otherbody++) {
        bforce(thisbody, otherbody, deltaf); // This function computes the pairwise force of otherbody on thisbody
        fx[thisbody] += deltaf[0]; // Add x component of force to thisbody
        fy[thisbody] += deltaf[1]; // Add y component of force to thisbody
        fz[thisbody] += deltaf[2]; // Add z component of force to thisbody
      }
    }

    // Now move the bodies (assumes constant acceleration during the timestep)
    for (thisbody=0; thisbody<n; thisbody++) {
      ax = fx[thisbody]/mass[thisbody]; // Compute x-direction acceleration of thisbody
      ay = fy[thisbody]/mass[thisbody]; // Compute y-direction acceleration of thisbody
      az = fz[thisbody]/mass[thisbody]; // Compute z-direction acceleration of thisbody
      vavgx = vx[thisbody] + dt2*ax; // Compute average x velocity of thisbody
      vavgy = vy[thisbody] + dt2*ay; // Compute average y velocity of thisbody
      vavgz = vz[thisbody] + dt2*az; // Compute average z velocity of thisbody
      x[thisbody] = x[thisbody] + dt*vavgx; // Compute new x position of thisbody
      y[thisbody] = y[thisbody] + dt*vavgy; // Compute new y position of thisbody
      z[thisbody] = z[thisbody] + dt*vavgz; // Compute new z position of thisbody
      vx[thisbody] += dt*ax; // Compute x velocity of thisbody at end of timestep
      vy[thisbody] += dt*ay; // Compute y velocity of thisbody at end of timestep
      vz[thisbody] += dt*az; // Compute z velocity of thisbody at end of timestep   
   }
  }

  // gather results
  // gather numbers
  MPI_Gather(&n, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  recv_offsets[0] = 0;
  for(i = 1; i < size; i++){
    recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
  }
  // gatherv x
  MPI_Gatherv(x, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&x, &buf);
  
  // gatherv y
  MPI_Gatherv(y, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&y, &buf);
  
  // gatherv z
  MPI_Gatherv(z, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&z, &buf);
  
  // gatherv mass
  MPI_Gatherv(mass, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&mass, &buf);

  // gatherv vx
  MPI_Gatherv(vx, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&vx, &buf);
  
  // gatherv vy
  MPI_Gatherv(vy, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&vy, &buf);
  
  // gatherv vz
  MPI_Gatherv(vz, n, MPI_DOUBLE, buf, recv_counts, recv_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  swap(&vz, &buf);



  if(rank == 0){
    n = N;
    output(K); // Final output
    timing(&etime1, &cptime);
    etime = etime1 - etime0;
    printf ("\nTime for %d timesteps with %d bodies: %9.4f seconds\n", K, N, etime);
  }


  // Free arrays
  free(mass);
  free(x);
  free(y);
  free(z);
  free(vx);
  free(vy);
  free(vz);
  free(fx);
  free(fy);
  free(fz);

  free(bx);
  free(by);
  free(bz);
  free(bmass);
  free(buf);  
  MPI_Finalize();
}
