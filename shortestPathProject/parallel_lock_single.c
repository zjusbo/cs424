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
  struct s_adj_node* next;
} adj_node;


/*Moore's Algorithm*/
int N; // number of vertices
int* sources; // source node
int num_sources;

omp_lock_t dist_lock, queue_lock;
adj_node** adj_listhead;

void process(int vi, int* dist, queue* pq, bool* isInQueue){
  adj_node* vj_p = adj_listhead[vi];
  
  // Loop over edges out of vi
  while(vj_p){
    int vj = vj_p -> vertex;
    int newdist_vj;
    bool needToUpdate;
      
    // Do this if new distance is smaller
    { 
      newdist_vj = dist[vi] + vj_p->weight;// new distance throught vi
      needToUpdate = (newdist_vj < dist[vj]) || (dist[vj] == INF); 
    }
    {
    if(needToUpdate){
      
      omp_set_lock(&queue_lock);
      dist[vj] = newdist_vj;
      if(isInQueue[vj] == false){
        q_enqueue(vj, pq);
        isInQueue[vj] = true;
      }
      omp_unset_lock(&queue_lock);
    }
    vj_p = vj_p -> next;
    }
  }
}

int moore(int source){
  // distance between source vertex and current vertex
  int* dist;
  bool* isInQueue;
  queue q;
  // Initialize
  dist =(int *) malloc((N+1) * sizeof(int));     
  isInQueue =(bool *) malloc((N+1) * sizeof(bool));     
  for(int i = 1; i <= N; i++) dist[i] = INF;
  for(int i = 1; i <= N; i++) isInQueue[i] = false; // queue is empty

  q_init(&q);
  dist[source] = 0;
  q_enqueue(source, &q); 
  // Loop over entries in queue
  while(true){
  #pragma omp parallel shared(dist, adj_listhead, q, dist_lock, queue_lock) 
  #pragma omp single
  while(!q_isEmpty(&q)){
    int vi;
    omp_set_lock(&queue_lock);
    vi = q_dequeue(&q);
    isInQueue[vi] = false;
    omp_unset_lock(&queue_lock);
    #pragma omp task 
    {
      process(vi, dist, &q, isInQueue);
    }
  } // All done
  // implicit barrier
  // all tasks should be finished below this line
  if(q_isEmpty(&q))
    break;
  }
  if(DEBUG){
    printf("source = %d, ", source);
    printf("%d %d %d", dist[1], dist[N-1], dist[N]);
    printf("\n");
  }
  free(dist);
  free(isInQueue);
}

void adj_list_add(int vertex, int weight, adj_node** adj_head){
   adj_node* node = (adj_node*)malloc(sizeof(adj_node));
   node -> vertex = vertex;
   node -> weight = weight; 
   node -> next = NULL;
   while(*adj_head != NULL){
     adj_head = &((*adj_head) -> next);
   }
   *adj_head = node;
}

void print_adj_list(adj_node** adj_head, int n){
  for(int i = 1; i <= n; i++){
    printf("%d: ", i);
    adj_node * node = adj_head[i];
    while(node != NULL){
      printf("%d -> ", node->vertex);
      node = node -> next;
    }
    printf("NULL\n");
  }
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
     adj_listhead = (adj_node**)malloc(sizeof(adj_node**) * (N + 1));
     for(int i = 1; i <= N; i++) adj_listhead[i] = NULL;
   }else if(line[0] == 'a'){
     // data line  
     sscanf(line, "%s %d %d %d", tmp, &start, &end, &weight);
     adj_list_add(end, weight, &adj_listhead[start]);
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
  omp_init_lock(&dist_lock);
  omp_init_lock(&queue_lock);

  for(int i = 0; i < num_sources; i++){
  //  if(DEBUG)
//      printf("Computing source %d\n", sources[i]);
    moore(sources[i]);
  }
  for(int i = 1; i <=N; i++){
    adj_node* node = adj_listhead[i];
    while(node){
      adj_node* next = node->next;
      free(node);
      node = next;
    }
  }
  free(sources);
  timing(&wct1, &cput); //get the end time
  total_time = wct1 - wct0;
  printf("Message printed by master: Total elapsed time is %f seconds.\n",total_time);
}
