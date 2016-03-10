#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"

/*
  This is the serial program for CPSC424/524 Assignment #4.
  
  Author: Andrew Sherman, Yale University

  Date: 2/23/2013

*/

// Globals so that we don't need to pass them to functions
double dt, dt2;
int N, K;
double *mass, *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz;

// Function to compute the forces between a pair of bodies
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
  for (thisbody=0; thisbody<N; thisbody++) {
    cmassx += mass[thisbody]*x[thisbody];
    cmassy += mass[thisbody]*y[thisbody];
    cmassz += mass[thisbody]*z[thisbody];
    tmass  += mass[thisbody];
    tvelx += vx[thisbody];
    tvely += vy[thisbody];
    tvelz += vz[thisbody];
  }
  printf("\n     Center of Mass:   (%e, %e, %e)\n", cmassx/tmass, cmassy/tmass, cmassz/tmass);
  printf(  "     Average Velocity: (%e, %e, %e)\n", tvelx/N, tvely/N, tvelz/N);
}


// Main function
int main(int argc, char **argv) {

  int i, ts, thisbody, otherbody;
  double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
  // Serial Timing Functions: See the timer.h include file in this directory
  double etime, etime0, etime1, cptime;


  // Read basic inputs   
  scanf("%d\n",&N);
  scanf("%d\n",&K);
  scanf("%le\n",&dt);
  
  dt2 = dt/2.;

  // Allocate arrays
  mass = calloc(N, sizeof(double)); // mass[i] is mass of body i
  x = calloc(N, sizeof(double)); // x[i] is x position of body i
  y = calloc(N, sizeof(double)); // y[i] is y position of body i
  z = calloc(N, sizeof(double)); // y[i] is z position of body i
  vx = calloc(N, sizeof(double)); // vx[i] is x velocity of body i
  vy = calloc(N, sizeof(double)); // vy[i] is y velocity of body i
  vz = calloc(N, sizeof(double)); // vz[i] is z velocity of body i
  fx = calloc(N, sizeof(double)); // fx[i] is x force on body i
  fy = calloc(N, sizeof(double)); // fy[i] is y force on body i
  fz = calloc(N, sizeof(double)); // fz[i] is z force on body i

  // Read initial conditions
  for (i=0; i<N; i++) scanf("%le\n",&mass[i]);
  for (i=0; i<N; i++) scanf("%le %le %le\n",&x[i],&y[i],&z[i]);
  for (i=0; i<N; i++) scanf("%le %le %le\n",&vx[i],&vy[i],&vz[i]);

  // Start Timer
  timing(&etime0,&cptime);

  // Timestep Loop
  for (ts=0; ts<K; ts++) {
    if (ts%128 == 0) output(ts); // Print output if necessary

    // Initialize forces on bodies
    for (thisbody=0; thisbody<N; thisbody++) {
      fx[thisbody] = 0.;
      fy[thisbody] = 0.;
      fz[thisbody] = 0.;
    }

    // Compute all pairwise interbody forces (Note that we take advantage of symmetry.)
    for (thisbody=0; thisbody<N; thisbody++) {
      for (otherbody=thisbody+1; otherbody<N; otherbody++) {
	force(thisbody, otherbody, deltaf); // This function computes the pairwise force of otherbody on thisbody
	fx[thisbody] += deltaf[0]; // Add x component of force to thisbody
	fy[thisbody] += deltaf[1]; // Add y component of force to thisbody
	fz[thisbody] += deltaf[2]; // Add z component of force to thisbody
	fx[otherbody] -= deltaf[0]; // Subtract x component of force from otherbody
	fy[otherbody] -= deltaf[1]; // Subtract y component of force from otherbody
	fz[otherbody] -= deltaf[2]; // Subtract z component of force from otherbody
      }
    }

    // Now move the bodies (assumes constant acceleration during the timestep)
    for (thisbody=0; thisbody<N; thisbody++) {
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
  output(K); // Final output

  timing(&etime1, &cptime);
  etime = etime1 - etime0;
  printf ("\nTime for %d timesteps with %d bodies: %9.4f seconds\n", K, N, etime);

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
}

