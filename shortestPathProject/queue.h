#ifndef _QUEUE
#define _QUEUE
#define _GNU_SOURCE // supress getline warning
#include <stdio.h>
#include <stdbool.h>


typedef struct s_q_node{
  int elem;
  struct s_q_node * next;
}q_node;

typedef struct s_queue{
  q_node * head, * tail;
}queue;

void q_init(queue * q);
bool q_isEmpty(queue* q);
void q_enqueue(int elem, queue * q);
int q_dequeue(queue * q);

#endif
