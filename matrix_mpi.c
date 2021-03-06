#include <stdio.h>
#include <stdlib.h>
#include "gen_matrix.h"
#include "my_malloc.h"
#include "mpi.h"


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

main(int argc, char *argv[]) {
  double **r;
  double **result;
  int i;
  int num_arg_matrices;

  int rank, size;
  MPI_Status status;
    
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  double start, end;
  start = MPI_Wtime();

  if (argc != 4) {
    printf("usage: debug_perf test_set matrix_dimension_size\n");
    exit(1);
  }
  int debug_perf = atoi(argv[1]);
  int test_set = atoi(argv[2]);
  matrix_dimension_size = atoi(argv[3]);
  num_arg_matrices = init_gen_sub_matrix(test_set);
  
  // allocate arrays
  r = (double **)my_malloc(sizeof(double *) * num_arg_matrices);
  result = (double **)my_malloc(sizeof(double *) * 2);
  result[0] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
  result[1] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);

  // get sub matrices
  for (i = 0; i < num_arg_matrices; ++i) {
    r[i] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
    if (gen_sub_matrix(0, test_set, i, r[i], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
      printf("inconsistency in gen_sub_matrix\n");
      exit(1);
    }
  }  

  // start MPI implementation
  
  
  // perform matrix multiplies of r[0] and r[1]
   MPI_Barrier(MPI_COMM_WORLD);
   int x,y,k;
   for(x = rank; x < matrix_dimension_size; x = x +size) //divide the task in multiple processes
    {
        for(y = 0; y < matrix_dimension_size; y++)
        {
            double sum=0;
            for(k = 0; k < matrix_dimension_size; k++)
            {
                sum = sum + r[0][x*matrix_dimension_size+ k] * r[1][k*matrix_dimension_size+y];
            }
            result[0][x*matrix_dimension_size + y] = sum;
        }
    }
    
    if(rank != 0)
    {
        for(x = rank; x < matrix_dimension_size; x = x+size)
        MPI_Send(&result[0][x*matrix_dimension_size + 0], matrix_dimension_size, MPI_DOUBLE, 0, 10+x, MPI_COMM_WORLD);//send calculated rows to process with rank 0
    }
    
    if(rank == 0)
    {
        for(y = 1; y < size; y++)
        {
            for(x = y; x < matrix_dimension_size; x = x+size)
            {
                MPI_Recv(&result[0][x*matrix_dimension_size+0], matrix_dimension_size, MPI_DOUBLE, y, 10+x, MPI_COMM_WORLD, &status);//receive calculated rows from respective process
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*
    int n = 0;
    if (debug_perf == 0) {
      // print each of the sub matrices
     for (i = 0; i < num_arg_matrices; ++i) {
        printf("argument matrix %d\n", i);
       print_matrix(r[i], matrix_dimension_size);
      }
      printf("result matrix\n");
      print_matrix(result[n], matrix_dimension_size);
    } 
    */
    //MPI_Bcast(&result[0][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE,0,MPI_COMM_WORLD);
    int n =0;
    i=2;
    while( i < num_arg_matrices) {
      MPI_Bcast(&result[n][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      int x,y,k;
      for(x = rank; x < matrix_dimension_size; x = x+size) //divide the task in multiple processes
      {
          
          for(y = 0; y < matrix_dimension_size; y++)
          {
              double sum=0;
              for(k = 0; k < matrix_dimension_size; k++)
              {
                 sum = sum + result[n][x*matrix_dimension_size+ k] * r[i][k*matrix_dimension_size+y];
              }
              result[n ^ 0x1][x*matrix_dimension_size + y] = sum;
          }
      }
    
      if(rank != 0)
      {
          for(x = rank; x < matrix_dimension_size; x = x+size)
          MPI_Send(&result[n ^ 0x1][x*matrix_dimension_size + 0], matrix_dimension_size, MPI_DOUBLE, 0, 10+x, MPI_COMM_WORLD);//send calculated rows to process with rank 0
      }
      
      if(rank == 0)
      {
         for(y = 1; y < size; y++)
          {
              for(x = y; x < matrix_dimension_size; x = x+size)
              {
                  MPI_Recv(&result[n ^ 0x1][x*matrix_dimension_size+0], matrix_dimension_size, MPI_DOUBLE, y, 10+x, MPI_COMM_WORLD, &status);//receive calculated rows from respective process
             }
          }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      n = n ^ 0x1;
      i++;
      
    }
 
  if(rank ==0 ){
    end = MPI_Wtime();
    if (debug_perf == 0) {
      // print each of the sub matrices
      for (i = 0; i < num_arg_matrices; ++i) {
        printf("argument matrix %d\n", i);
        print_matrix(r[i], matrix_dimension_size);
      }
      printf("result matrix\n");
      print_matrix(result[n], matrix_dimension_size);
     // printf("\n");
     // printf("Runtime = %f\n", end-start);
    } else {
      double sum = 0.0;

      for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
        sum += result[n][i];
      }
      printf("%f\n", sum);
    }
  }
  MPI_Finalize();
}

