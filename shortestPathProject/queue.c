#include "queue.h"
#include <stdlib.h>

bool q_isEmpty(queue* q){
  return (q -> head == NULL) && (q -> tail == NULL);
}

void q_init(queue * q){
  q->head = q->tail = NULL;
}

// TODO check if elem already exists in queue?
void q_enqueue(int elem, queue * q){
  q_node* node = (q_node *)malloc(sizeof(q_node));
  node -> elem = elem;
  node -> next = NULL;
  if(q_isEmpty(q)){
    q -> head = q -> tail = node;
  }else{
    q -> tail -> next = node;
    q -> tail = node;
  }
}

int q_dequeue(queue * q){
  int elem;
  if(q_isEmpty(q)){
    printf("[error] Queue is empty\n");
    return -1; // ERROR , queue is empty
  }
  q_node* node = q -> head;
  q -> head = q -> head -> next;
  if(q -> head  == NULL) q -> tail = NULL;
  elem = node->elem;
  free(node);
  return elem;
}
