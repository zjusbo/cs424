#include <stdio.h>
#include <stdlib.h>
int cal_block_size(int, int, int);
int main(int argc, char ** argv){
	int i;
	int N = 100;
	int num_nodes = 10;
	for(i = 0; i < num_nodes; i++){
		int size = cal_block_size(N, i, num_nodes);
		printf("%d, %d\n", i, size);
	}

}

int * _block_size = NULL;
int cal_block_size(int N, int rank, int num_nodes){
  if(_block_size != NULL){
    // return previously calculated value directly
    return _block_size[rank];
  }
  else{
    // calculate block size for all processes
    int total_len = (1 + N) * N / 2;
    int average_len = total_len / num_nodes;
    _block_size = (int* )malloc(sizeof(int) * num_nodes);
    int node_idx = 0;
    int i;
    int sum = 0;
    int row_idx = 0;
    for(i = 1; i <= N; i++){
      sum += i;
      if(sum > average_len){ // not the last node
        i--;
        _block_size[node_idx++] = i - row_idx;
        row_idx = i;
        sum = 0;
      }
      if(node_idx == num_nodes - 1){
        // last node
        _block_size[node_idx] = N - row_idx;
        break;
      }
    }
  }
}