#include <stdio.h>
#include <stdlib.h>
int cal_block_size(int, int, int);
int main(int argc, char ** argv){
	int i;
	int N = atoi(argv[1]);
	int num_nodes = atoi(argv[2]);
    int sum = 0;
	int size = cal_block_size(N, i, num_nodes);
	return 0;
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
        printf("%dth node, height: %d, len: %d\n", node_idx, i - row_idx, sum - i);
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
    return cal_block_size(N, rank, num_nodes);
  }
}
