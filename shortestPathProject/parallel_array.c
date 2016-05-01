#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <omp.h>
#include "queue.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#define DEBUG 0
#define INF -1

typedef struct s_adj_node{
  int vertex;
  int weight;
} adj_node;



/*Moore's Algorithm*/
int N; // number of vertices
int* sources; // source node
int num_sources;

adj_node** adj_listhead;
int* adj_list_length;
void process(int vi, int* dist, queue* pq){
  adj_node* vj_p_base = adj_listhead[vi];
  // Loop over edges out of vi
  #pragma omp parallel
  #pragma omp for schedule(dynamic)
  for(int i = 0; i < adj_list_length[vi]; i++){
  //  printf("Thread %d: i = %d, vj_p -> vertex = %d\n", omp_get_thread_num(), i, vj_p->vertex);
    int vj = vj_p_base[i].vertex;
    int weight = vj_p_base[i].weight;
  //  printf("Thread %d: i = %d, vj = %d\n", omp_get_thread_num(), i, vj); 
    int newdist_vj;
    bool needToUpdate;
 
    // Do this if new distance is smaller
    newdist_vj = dist[vi] + weight;// new distance throught vi
    needToUpdate = (newdist_vj < dist[vj]) || (dist[vj] == INF); 

    if(needToUpdate){
      dist[vj] = newdist_vj;
      #pragma omp critical(queue)
      q_enqueue(vj, pq);
    }
   // printf("Thread %d: i = %i, end of loop\n", omp_get_thread_num(), i);
  }
}

int moore(int source){
  // distance between source vertex and current vertex
  int* dist;
  queue q;
  // Initialize
  dist =(int *) malloc((N+1) * sizeof(int));     
  for(int i = 1; i <= N; i++) dist[i] = INF;
  q_init(&q);
  dist[source] = 0;
  q_enqueue(source, &q); 
  // Loop over entries in queue
  while(!q_isEmpty(&q)){
    int vi;
    vi = q_dequeue(&q);
    process(vi, dist, &q);
  } // All done
  // implicit barrier
  // all tasks should be finished below this line
  if(DEBUG){
    printf("source = %d, ", source);
    printf("%d %d %d", dist[1], dist[N-1], dist[N]);
    printf("\n");
  }
  free(dist);
}



void readGraph(char* filename){
 FILE* fp;
 int num_nodes, num_edges;
 int start, end, weight;
 char tmp[10];
 char* line;
 size_t len = 100;
 fp = fopen(filename, "r");
 if(fp == NULL){
   printf("can not open graph file %s\n", filename);
   return ;
 }
 line = (char*)malloc(sizeof(char) * len);
 while(getline(&line, &len, fp)	!= -1){
   if(line[0] == 'c') continue; // comment line
   if(line[0] == 'p'){
     // program line
     sscanf(line, "%s %s %d %d", tmp, tmp, &num_nodes, &num_edges);
     N = num_nodes;
     // initialization
     adj_listhead = (adj_node**)malloc(sizeof(adj_node*) * (N + 1));
     adj_list_length = (int*)malloc(sizeof(int) * (N + 1));
     for(int i = 1; i <= N; i++) adj_listhead[i] = (adj_node*)malloc(sizeof(adj_node) * N);
     for(int i = 1; i <= N; i++) adj_list_length[i] = 0;
   }else if(line[0] == 'a'){
     // data line  
     sscanf(line, "%s %d %d %d", tmp, &start, &end, &weight);
     int idx = adj_list_length[start];
     adj_listhead[start][idx].vertex = end;
     adj_listhead[start][idx].weight = weight;
     //adj_list_add(end, weight, &adj_listhead[start]);
     adj_list_length[start]++;
   }
 }
 fclose(fp);
 free(line);
}

bool readSource(char* filename){
  FILE* fp;
  char tmp[10];
  char* line;
  size_t len = 100; 
  int source;
  int idx;
  fp = fopen(filename, "r"); 
  if(fp == NULL){
    printf("can not open source file %s\n", filename);
    return false;
  }
  line = (char*)malloc(sizeof(char) * len);
  while(getline(&line, &len, fp) != -1){
    if(line[0] == 'c') continue;
    else if(line[0] == 'p'){
      sscanf(line, "%s %s %s %s %d", tmp, tmp, tmp, tmp, &num_sources);
      sources = (int*)malloc(sizeof(int) * num_sources);
      idx = 0;
    }else if(line[0] == 's'){
      sscanf(line, "%s %d", tmp, &source);
      sources[idx++] = source;
    }
  }
  fclose(fp);
  free(line);
}

int main(int argc, char **argv ) {

  /*
    This is the shortest path project for CPSC424/524.

    Author: Bo Song, Yale University

    Date: 4/25/2016

    Credits: This program is based on the description provided by Andrew Sherman
  */

  double wct0, wct1, total_time, cput;
  char* sourceFile, * graphFile;
  #pragma omp parallel
  printf("num of threads = %d\n", omp_get_num_threads());
  if(argc != 3){
    printf("serial <graphfile> <sourcefile>\n");
    return -1;
  }
  graphFile = argv[1];
  sourceFile = argv[2];
  timing(&wct0, &cput);  
  printf("reading graph...\n");
  readGraph(graphFile);
  printf("reading source...\n");
  readSource(sourceFile);
  // print_adj_list(adj_listhead, N);
  for(int i = 0; i < num_sources; i++){
  //  if(DEBUG)
  //    printf("Computing source %d\n", sources[i]);
    moore(sources[i]);
  }
  timing(&wct1, &cput); //get the end time
  total_time = wct1 - wct0;
  printf("Message printed by master: Total elapsed time is %f seconds.\n",total_time);
}
