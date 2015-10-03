#include <stdio.h>
#include <stdlib.h>
#include "gen_matrix.h"
#include <cilk/cilk.h>
#include "my_malloc.h"
#include <time.h>

void copyMatrix(double* result, double * a, int size){
    int i;
    for (i = 0 ; i < size * size ; i++){
       result[i] = a[i];
    }
}

void multMatrices( double* A, double* B, double* C, int size, int length, int ctr)
{
  int j, h, k ;
  if(length == 1) 
  {

    for(j=0; j<size; j++)
    {
      k = ctr;
      C[ctr*size+j] = 0;
      for(h = 0; h < size; h++)
      {
        C[ctr*size+j] += A[k*size+h] * B[h*size+j];
      }   
    }

  }
  else
  {
    cilk_spawn multMatrices(A, B, C, size, length/2, ctr);
    cilk_spawn multMatrices(A, B, C, size,length/2, ctr + length/2);
    if(length % 2 != 0)
    {
      cilk_spawn multMatrices(A, B, C, size, 1, ctr + length - 1);
    }
    cilk_sync;
  }
}

void mm(double * C, double * A, double *B, int n){
  int  split, i, ctr = 0, nchunk, chunk;
  for(split= n/2+1; (split <= n ) && (n % split != 0); split++){
  }
  
  chunk = n/split;
  
  for(i=0; i<split; i++)
  {
    if(ctr == 0)
    {
      cilk_spawn multMatrices(A, B, C, n, chunk, ctr);
      cilk_sync;
      ctr += chunk;
    }
    else if(i == split - 1)
    {
      nchunk = n - chunk * i;
      ctr = chunk * i;
      cilk_spawn multMatrices(A, B, C, n, nchunk, ctr);  
      cilk_sync;
    }
    else
    {
      cilk_spawn multMatrices(A, B, C, n, chunk, ctr); 
      cilk_sync;
      ctr += chunk;
    }
  }
}

void print_matrix(double *result, int dim_size) {
  int x, y;
  for (y = 0; y < dim_size; ++y) {
    for (x = 0; x < dim_size; ++x) {
      printf("%f ", result[y * dim_size + x]);
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char *argv[]) {
 // clock_t t = clock();
  double **r;
  double **result;
  double **back_up;
  int  i;
  int num_arg_matrices;
  int debug_perf, test_set, tmp , gap , qq, k;
  double sum ;

  // test running time
  double time_exc ;
  clock_t t;
  t = clock();
  if (argc != 4) {
    printf("usage: debug_perf test_set matrix_dimension_size\n");
    exit(1);
  }
  debug_perf = atoi(argv[1]);
  test_set = atoi(argv[2]);
  matrix_dimension_size = atoi(argv[3]);
  num_arg_matrices = init_gen_sub_matrix(test_set);
  
  // allocate arrays
  r = (double **)my_malloc(sizeof(double *) * num_arg_matrices);
  result = (double **)my_malloc(sizeof(double *) * num_arg_matrices/2);
  back_up = (double **)my_malloc(sizeof(double *) * num_arg_matrices/2);

  if(num_arg_matrices ==1){
         result[0] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
          if (gen_sub_matrix(0, test_set, 0, result[0], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
              printf("inconsistency in gen_sub_matrix\n");
              exit(1);
          }
          if (debug_perf == 0) {
            // print each of the sub matrices
            printf("argument matrix %d\n", i);
            print_matrix(result[0], matrix_dimension_size);
          }
    
         printf("result matrix\n");
         print_matrix(result[0], matrix_dimension_size);

         } else {
            sum = 0.0;

            for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
              sum += result[0][i];
            }
            printf("%f\n", sum);
          }
    my_free(result[0]);
    return 0;
  }

  for(qq  = 0 ; qq < num_arg_matrices/2; qq ++){
       result[qq] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
       back_up[qq] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
  }
  // get sub matrices
  for (i = 0; i < num_arg_matrices; ++i) {
    r[i] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
    if (gen_sub_matrix(0, test_set, i, r[i], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
      printf("inconsistency in gen_sub_matrix\n");
      exit(1);
    }
  }  

  
  // perform matrix multiplies
  for(k=0;k<num_arg_matrices/2;k++){
          cilk_spawn mm(result[k], r[k], r[k+1], matrix_dimension_size);
  }
  cilk_sync;
  
  if(num_arg_matrices%2){
     copyMatrix(back_up[num_arg_matrices/2-1], result[num_arg_matrices/2-1],matrix_dimension_size);
     cilk_spawn mm(result[num_arg_matrices/2-1], back_up[num_arg_matrices/2-1], r[num_arg_matrices-1], matrix_dimension_size);
     cilk_sync;
   }
  // calculate from buttom up
  
  tmp = num_arg_matrices/2;
  
  gap = 1 ;
  while(tmp!=1){
    int p;
    for(p=0;p<tmp/2;p++){
          copyMatrix(back_up[p], result[p], matrix_dimension_size);
          cilk_spawn mm(result[p], back_up[p], result[p+gap] , matrix_dimension_size);
    }
    cilk_sync;
    if(tmp%2){
        copyMatrix(back_up[tmp/2-1], result[tmp/2-1], matrix_dimension_size);
        cilk_spawn mm(result[tmp/2-1], back_up[tmp/2-1], result[tmp-1] , matrix_dimension_size);
        cilk_sync;
    }
    tmp = tmp/2;
    gap = gap*2;
  }

  t = clock()-t;
  if (debug_perf == 0) {
    // print each of the sub matrices
    
    for (i = 0; i < num_arg_matrices; ++i) {
      printf("argument matrix %d\n", i);
      print_matrix(r[i], matrix_dimension_size);
    }
    
    printf("result matrix\n");
    print_matrix(result[0], matrix_dimension_size);
   
    //time_exc =  ((double) t ) / CLOCKS_PER_SEC ;

   // printf("running time of cilk %f \n" , time_exc);
  } else {
    sum = 0.0;

    for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
      sum += result[0][i];
    }
    printf("%f\n", sum);
  }

  int ff =0;
  for(ff= 0 ; ff < num_arg_matrices ;mm++){
      my_free(r[ff]);
      if(ff%2 == 0){
        my_free(result[ff]);
        my_free(back_up[ff]);
      }
  }
}

