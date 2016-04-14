//
// Sequential Mandelbrot Calculation
//

/* 
   This program is based on a program provided in the solutions guide 
   for the text "Parallel Programming" by Barry Wilkinson and Michael Allen.
   The X Windows code and the basic computation in worker() have been retained,
   but the entire program has been restructured and the implementation of the
   dynamic task pool has been completely rewritten. In addition, the algorithm
   was modified to allow for randomization of the order of the computation.
   Finally command line arguments were added to vary the number of tasks, 
   make the randomization optional, and to permit refocusing of the
   computation onto different regions of the complex plane.

   Usage: mand [ntasks [irand [swr swi ner nei] ] ]
          ntasks: Ignored for sequential version (ntasks always 1)

          irand  == 0 or omitted: natural row order (default)
	  irand  != 0: randomized row order

          (swr,swi) and (ner, nei) delimit lower left and upper right corners
	  of focus region in the complex plane. If omitted, the focus region
	  has corners (-1.5,-1.) and (1., 1.).

   Programmer: Andrew Sherman, Yale University
   Date: February 2010
*/

#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#define X_RESN 1300// x resolution
#define Y_RESN 800// y resolution

#define MAXTASKS 800// Max number of tasks
#define MAXITNS 10000// Max number of Mandelbrot iterations per point

typedef struct complextype {
  float real, imag;
} Compl;

char USAGE[]=
  "Usage: mand [ntasks [irand [swr swi ner nei] ] ]\n\
          ntasks: Ignored for sequential version (ntasks always 1)\n\n\
          irand  ==0 or omitted: natural row order (default)\n\
          irand  !=0: randomized row order\n\n\
          (swr,swi) and (ner, nei) delimit lower left and upper right corners\n\
	  of focus region in the complex plane. If omitted, the focus region\n\
	  has corners (-1.5,-.5) and (1., 1.).";

int k;
double time1, time2, total_time, cptime;

float deltax, deltay;
float sw[2]={-1.5, -1.};
float ne[2]={1., 1.};

int rows[Y_RESN]; // row ordering for tasks

int irand=0;

int i, j, numitns;
Compl z, c;
float lengthsq, temp;

/* X Windows Variables */
GC gc;
Display *display;
Window win;


int main( int argc, char **argv ){

  setup(argc,argv);
  mandelbrot();

  return 0;
}

int mandelbrot() {

  /* START PROBLEM SOLUTION HERE */

  timing(&time1, &cptime);
      
  for(j = 0; j < Y_RESN; j++) {
    //work on pixel row in rows[j]
    c.imag = sw[1] + (float) (Y_RESN - rows[j]) * deltay;
    for(i=0; i < X_RESN; i++){
      
      c.real = sw[0] + (float) i * deltax;
      
      z.real = z.imag = 0.0;
      numitns = 0;
      do { /* iterate for pixel color */
	temp = z.real*z.real - z.imag*z.imag + c.real;
	z.imag = 2.0*z.real*z.imag + c.imag;
	z.real = temp;
	lengthsq = z.real*z.real+z.imag*z.imag;
	numitns++;
      } while ((lengthsq < 4.0) && (numitns < MAXITNS));
      XSetForeground (display, gc, (100*(MAXITNS-numitns)));
      XDrawPoint (display, win, gc, i, rows[j]); // Draw the pixel
    }
  }
  
  timing(&time2, &cptime);
  total_time = time2 - time1;
  
  XFlush (display);
  
  sleep (10);
  
  printf("Total Time: %f seconds\n",total_time);
  
}

int setup ( int argc, char **argv ){
int c;
  /* Check Command Line Arguments */

  while ((c = getopt (argc,argv,"+h")) != -1)
    switch (c)
      {
      case 'h':
        printf("%s\n", USAGE);
        exit(0);
        break;
      }

  if (argc>=2) {
    if (argc>=3) {
      irand = atoi(argv[2]);
      if (argc==7) {
	sw[0] = atof(argv[3]);
	sw[1] = atof(argv[4]);
	ne[0] = atof(argv[5]);
	ne[1] = atof(argv[6]);
      }
      else if (argc>3) {
	printf("Incorrect number of arguments.\n%s\n", USAGE);
	exit (-1);
      }
    }
  }
  
  // Calculate the distance in complex plane between pixels  
  deltax = (ne[0]-sw[0]) / (float) X_RESN;
  deltay = (ne[1]-sw[1]) / (float) Y_RESN;

  // Initialize pixel row ordering 
  //    (tasks are groups of pixel lines contained in consecutive entries of rows[])
  
  // Natural: rows ordered from 0 to Y_RESN-1.
  for(i=0; i<Y_RESN; i++) rows[i] = i;
  
  // Randomized: Randomly permute the order of rows
  if (irand) {
    srand(100*time(NULL));
    for(i=0; i<Y_RESN-1; i++){
      temp = (double)rand() / ((double)(RAND_MAX) + 1.);
      j = (int) (temp * (double)(Y_RESN-i)); 
      if (j > 0) {
	k = rows[i];
	rows[i] = rows[i+j];
	rows[i+j] = k;
      }
    }
  }

  setupX();
}


int setupX() {
/* X Windows Variables */
   unsigned int width, height,/* window size */
     x, y, /* window position */
     border_width,/* border width in pixels */
     display_width, display_height,/* size of screen */
     screen; /* which screen */
   char *window_name = "Mandelbrot Set", *display_name = NULL;
   unsigned long valuemask = 0;
   XGCValues values;
   XSizeHints size_hints;
   XSetWindowAttributes attr[1];
   
   /* Set up XWindows */
   if ( (display = XOpenDisplay (display_name)) == NULL ){
     fprintf (stderr, "drawon: cannot connect to X server %s\n",
	      XDisplayName (display_name) );
     exit (-1);
   }
   /* get screen size */
   screen = DefaultScreen (display);
   display_width = DisplayWidth (display, screen);
   display_height = DisplayHeight (display, screen);
   /* set window size */
   width = X_RESN;
   height = Y_RESN;
   /* set window position */
   x = 0;
   y = 0;
   /* create opaque window */
   border_width = 4;
   win = XCreateSimpleWindow (display, RootWindow (display, screen),
			      x, y, width, height, border_width,
			      BlackPixel (display, screen), WhitePixel (display, screen));
   size_hints.flags = USPosition|USSize;
   size_hints.x = x;
   size_hints.y = y;
   size_hints.width = width;
   size_hints.height = height;
   size_hints.min_width = 300;
   size_hints.min_height = 300;
   XSetNormalHints (display, win, &size_hints);
   XStoreName(display, win, window_name);
   /* create graphics context */
   gc = XCreateGC (display, win, valuemask, &values);
   XSetBackground (display, gc, WhitePixel (display, screen));
   XSetForeground (display, gc, BlackPixel (display, screen));
   XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
   attr[0].backing_store = Always;
   attr[0].backing_planes = 1;
   attr[0].backing_pixel = BlackPixel(display, screen);
   XChangeWindowAttributes(display, win, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);
   XMapWindow (display, win);
   XSync(display, 0);
   
 }
