//
// Parallel Mandelbrot with Dynamic Workload using non-blocking sends on workers 
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
          ntasks > 0: Dynamic task allocation using ntasks tasks
          ntasks == 0 or omitted: Static task allocation using nprocs-1 tasks (default)

          irand  == 0 or omitted: natural row order (default)
	  irand  != 0: randomized row order

          (swr,swi) and (ner, nei) delimit lower left and upper right corners
	  of focus region in the complex plane. If omitted, the focus region
	  has corners (-1.5,-1.) and (1., 1.).

   Programmer: Andrew Sherman, Yale University
   Date: February 2010
*/

#include "mpi.h"
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

#define MAXPROCS 64
#define MAXTASKS 800// Max number of tasks
#define MAXITNS 10000// Max number of Mandelbrot iterations per point

typedef struct complextype {
  float real, imag;
} Compl;

char USAGE[]=
  "Usage: mand [ntasks [irand [swr swi ner nei] ] ]\n\
          ntasks > 0: Dynamic task allocation using ntasks tasks\n\
          ntasks ==0 or omitted: Static task allocation using nprocs-1 tasks (default)\n\n\
          irand  ==0 or omitted: natural row order (default)\n\
          irand  !=0: randomized row order\n\n\
          (swr,swi) and (ner, nei) delimit lower left and upper right corners\n\
	  of focus region in the complex plane. If omitted, the focus region\n\
	  has corners (-1.5,-1.) and (1., 1.).";

int nprocs=0, myid=0, id=0;
char host[14]={'\0'};
MPI_Status status;
MPI_Request request;

int *kbuf, bufsize, bufsizei, k;
double time1, time2, total_time, cptime;

float deltax, deltay;
float sw[2]={-1.5, -1.};
float ne[2]={1., 1.};

int rows[Y_RESN]; // row ordering for tasks
char *bcastbuf; 
int bcastsize, size, bufpos;

int tasknum, ntasks, tasks_assigned, tasks_done, irand=0;
int task[3]; // array with tasknum, first rows index of task, last rows index of task

int i, j, numitns;
Compl z, c;
float lengthsq, temp;

/* X Windows Variables */
GC gc;
Display *display;
Window win;


int main( int argc, char **argv ){
  int setupcode;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (myid == 0) {
    setupcode = setup(argc,argv);
    //printf("setupcode = %d\n",setupcode);
    MPI_Bcast(&setupcode,1,MPI_INT,0,MPI_COMM_WORLD);
    if (setupcode==0) master();
  }
  else { 
    MPI_Bcast(&setupcode,1,MPI_INT,0,MPI_COMM_WORLD);
    if (setupcode==0) worker();
  }

  MPI_Finalize();
  return(0);
}

int master() {

  // Arrays only on master
  int taskstart[MAXTASKS+1];// first rows entry the tasks (index into rows[])
  int taskct[MAXPROCS]={0}; // counts number of tasks by each process
  double times[MAXPROCS]={0.}; // holds "useful" elapsed time in each process
  char hosts[MAXPROCS][14];
  
  /* MASTER-SPECIFIC SETUP HERE */

  // Allocate task result buffer
  bufsizei = X_RESN * ((Y_RESN+ntasks-1)/ntasks);
  bufsize = bufsizei * sizeof(int);
  kbuf = (int*) malloc(bufsize);
  
  // Set up tasks as consecutive entries of rows[], which contains pixel row numbers (possibly randomized)
  taskstart[0] = 0;
  for (i=0; i<ntasks; i++) taskstart[i+1] = taskstart[i] + (Y_RESN/ntasks) + (i < (Y_RESN%ntasks));
  
  /*BROADCAST GENERAL WORKER DATA */

  // Send workers problem data they require
  MPI_Pack_size(4, MPI_FLOAT, MPI_COMM_WORLD, &bcastsize);
  MPI_Pack_size(Y_RESN+1, MPI_INT, MPI_COMM_WORLD, &size);
  bcastsize += size;
  bcastbuf = (char*) malloc(bcastsize);
  
  bufpos = 0;
  MPI_Pack(sw, 2, MPI_FLOAT, bcastbuf, bcastsize, &bufpos, MPI_COMM_WORLD);
  MPI_Pack(&deltax, 1, MPI_FLOAT, bcastbuf, bcastsize, &bufpos, MPI_COMM_WORLD);
  MPI_Pack(&deltay, 1, MPI_FLOAT, bcastbuf, bcastsize, &bufpos, MPI_COMM_WORLD);
  MPI_Pack(&ntasks, 1, MPI_INT, bcastbuf, bcastsize, &bufpos, MPI_COMM_WORLD);
  MPI_Pack(rows, Y_RESN, MPI_INT, bcastbuf, bcastsize, &bufpos, MPI_COMM_WORLD);
  MPI_Bcast(bcastbuf, bcastsize, MPI_PACKED, 0, MPI_COMM_WORLD);
  
  free(bcastbuf);
    
  /* START PROBLEM SOLUTION HERE */

  // Barrier here for timing purposes
  MPI_Barrier(MPI_COMM_WORLD);
  timing(&time1, &cptime);
  
  tasks_assigned = nprocs - 1;
  if (ntasks < tasks_assigned) tasks_assigned = ntasks;
  
  // Send out initial tasks to every worker
  for(i=0; i<tasks_assigned; i++) {
    task[0] = i;
    task[1] = taskstart[i];
    task[2] = taskstart[i+1] - 1;
    MPI_Send(task, 3, MPI_INT, i+1, 1, MPI_COMM_WORLD);
    taskct[i+1]++;
  }
  // Shut down any unneeded workers
  for(id=tasks_assigned+1; id<nprocs; id++) {
    MPI_Send(task, 1, MPI_INT, id, 99, MPI_COMM_WORLD);
  }
  
  tasks_done = 0;
  while(tasks_done < ntasks){
    
    // Wait for any result
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    id = status.MPI_SOURCE; // Get the worker's id
    tasknum = status.MPI_TAG; // Task number came back in tag
    
    MPI_Irecv(kbuf, bufsizei, MPI_INT, id, tasknum, MPI_COMM_WORLD, &request);      
    // If there are unassigned tasks, send next one to this process
    if (tasks_assigned < ntasks) {
      task[0] = tasks_assigned;
      task[1] = taskstart[tasks_assigned];
      task[2] = taskstart[tasks_assigned+1] - 1;
      MPI_Send(task, 3, MPI_INT, id, 1, MPI_COMM_WORLD);
      taskct[id]++;
      tasks_assigned++;
    }
    // Otherwise, tell it to shut down
    else {
      MPI_Send(task, 3, MPI_INT, id, 99, MPI_COMM_WORLD);
    }
    // Wait for receive to complete
    MPI_Wait(&request, &status);
    
    k = 0;
    for(j=taskstart[tasknum]; j<=taskstart[tasknum+1]-1; j++) {
      for(i=0; i<X_RESN; i++, k++){ // kbuf[k] is 0 for points in the set
	XSetForeground (display, gc, (100*(MAXITNS-kbuf[k])));
	XDrawPoint (display, win, gc, i, rows[j]); // Draw the pixel
      }
    }
    tasks_done++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  timing(&time2, &cptime);
  total_time = time2 - time1;
  
  times[0] = total_time;
  for (id=1;id<nprocs;id++) {
    MPI_Recv(&times[id],1,MPI_DOUBLE,id,99,MPI_COMM_WORLD,&status);
    MPI_Recv(hosts[id],14,MPI_CHAR,id,99,MPI_COMM_WORLD,&status);
  }
  
  XFlush (display);
  
  sleep (10);
  
  printf("Total Time: %f seconds\n",total_time);
  printf("Task Counts & Times by Process:\n   proc host          tasks time (secs)\n   ---- ------------- ----- -----------\n");
  hosts[0][13] = '\0';
  gethostname(hosts[0], 13);
  for (id=0;id<nprocs;id++) printf("   %4d %s %4d  %f\n",id,hosts[id],taskct[id],times[id]);
  
return(0);
}

int worker() {
  int send_active=0, bufnum=0, task0;
  int *kbufs[2];
  MPI_Request requests[2];
  
  MPI_Pack_size(4, MPI_FLOAT, MPI_COMM_WORLD, &bcastsize);
  MPI_Pack_size(Y_RESN+1, MPI_INT, MPI_COMM_WORLD, &size);
  bcastsize += size;
  bcastbuf = (char*) malloc(bcastsize);
  
  MPI_Bcast(bcastbuf, bcastsize, MPI_PACKED, 0, MPI_COMM_WORLD); 
  
  bufpos = 0;
  MPI_Unpack(bcastbuf, bcastsize, &bufpos, sw, 2, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(bcastbuf, bcastsize, &bufpos, &deltax, 1, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(bcastbuf, bcastsize, &bufpos, &deltay, 1, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(bcastbuf, bcastsize, &bufpos, &ntasks, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(bcastbuf, bcastsize, &bufpos, rows, Y_RESN, MPI_INT, MPI_COMM_WORLD);
  
  free(bcastbuf);
  
  
  // Tasks are a whole number of horizontal pixel rows. Buffer is sized to hold a task of max size.
  bufsizei = X_RESN * ((Y_RESN+ntasks-1)/ntasks);
  bufsize = bufsizei*sizeof(int);
  kbufs[0] = (int*) malloc(bufsize);
  kbufs[1] = (int*) malloc(bufsize);
  
  // Barrier here for timing purposes
  MPI_Barrier(MPI_COMM_WORLD);
  timing(&time1, &cptime);
  // Loop forever as a worker
  while (1) {
    
    // At this point kbufs[bufnum] always points at an unused buffer

    // Get a task assignment
    MPI_Recv(task, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    if (status.MPI_TAG == 99) { // Tag=99 means we're all done
      if (send_active>1) {
	MPI_Wait(&requests[1-bufnum], &status); //Wait for prior send to complete
	send_active--;
      }
      break; 
    }

    k = 0;
    kbuf = kbufs[bufnum];
    for(j = task[1]; j <= task[2]; j++) {
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
	kbuf[k++] = numitns;
      }
    }

    task0 = task[0];
    MPI_Isend(kbufs[bufnum], bufsizei, MPI_INT, 0, task0, MPI_COMM_WORLD, &requests[bufnum]); // sending color
    send_active++;
    bufnum = 1 - bufnum;
    // We only have 2 buffers, so we wait if both are active.
    // bufnum will be the older one, so that's the one we wait for
    if (send_active>1) { 
      MPI_Wait(&requests[bufnum], &status); //Wait for prior send to complete
      send_active--;
    }
   
  }
  timing(&time2, &cptime);
  total_time = time2 - time1;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Send(&total_time, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
  gethostname(host, 13);
  MPI_Send(host,14,MPI_CHAR,0,99,MPI_COMM_WORLD);
  
  free(kbuf);
  
  return(0);
}

int setup ( int argc, char **argv ){
int c, setupcode;
  /* Check Command Line Arguments */

  while ((c = getopt (argc,argv,"+h")) != -1)
    switch (c)
      {
      case 'h':
        printf("%s\n", USAGE);
        return 1;
        break;
      }
  
  if (nprocs<2) {
    printf("Must use at least 2 MPI processes.\n");
    return -1;
  }
  if (argc>=2) {
    ntasks = atoi(argv[1]);
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
        return -1;
      }
    }
  }
  
  if (ntasks == 0) {
    ntasks = nprocs - 1;
    printf("Static task allocation using %d tasks.\n",ntasks);
  }
  else if (ntasks > 0) {
    printf("Dynamic task allocation using %d tasks.\n",ntasks);
  }
  else {
    printf("Invalid setting for ntasks.\n%s\n", USAGE);
    return -1;
  }
  if (ntasks>MAXTASKS) {
    printf("ntasks (=%d) is larger than the maximum permitted (=%d)\n",ntasks,MAXTASKS);
    return -1;
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

  setupcode = setupX();
  return setupcode;
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
     return -1;
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
   return (0);
 }
