#include <stdio.h>
#include <stdlib.h>
#include "gen_matrix.h"
#include "my_malloc.h"
#include "mpi.h"


void copyMatrix(double* result, double * a, int size){
    int i;
    for (i = 0 ; i < size * size ; i++){
       result[i] = a[i];
    }
}

// row major
void mm(double *result, double *a, double *b, int dim_size) {
  int x, y, k;
  for (y = 0; y < dim_size; ++y) {
    for (x = 0; x < dim_size; ++x) {
      double r = 0.0;
      for (k = 0; k < dim_size; ++k) {
	r += a[y * dim_size + k] *  b[k * dim_size + x];
      }
      result[y * dim_size + x] = r;
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

main(int argc, char *argv[]) {

  //=======================Common Data ================================================================//
  double **r;
  double **rp;
  double *result;
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
  
  //==================================Number of Processer is not enough=================================//
  if(size < num_arg_matrices/2){
           
           int avgmatrix = num_arg_matrices/size;
           int extraMatrix = num_arg_matrices - avgmatrix*size;
           
           //=============================Allocate Memory===============================================//
            r = (double **)my_malloc(sizeof(double *) * 2);
            r[0] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
            r[1] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
            result = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);


           // ==============================Calculate Result on Each Process=============================================//

            int dest;
            int offset = 0;
            int numOfMatrix = 0;
            for(dest = 0 ; dest< size; dest ++){
                numOfMatrix = (dest< extraMatrix) ? avgmatrix+1 : avgmatrix;
                if(dest==rank){
                     int tmp = numOfMatrix;
                     int i = 0;
                     if (gen_sub_matrix(0, test_set, offset+i, r[0], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                           printf("inconsistency in gen_sub_matrix\n");
                           exit(1);
                           }
                     if (gen_sub_matrix(0, test_set, offset+i+1, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                           printf("inconsistency in gen_sub_matrix\n");
                           exit(1);
                           }
                     mm(result,r[0],r[1],matrix_dimension_size);
                     for( i = 2 ; i <tmp ; i++){
                          copyMatrix(r[0],result,matrix_dimension_size);
                          if (gen_sub_matrix(0, test_set, offset+i, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                           printf("inconsistency in gen_sub_matrix\n");
                           exit(1);
                          }
                          mm(result,r[0],r[1],matrix_dimension_size);
                     }
                }
                offset+=numOfMatrix;
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
           //===============================Combine result from bottom up ========================================//
           int tmpsize = size;
           int previousNumOfProcess =1;
           int numOfprocessor = 2;
           int offsetRow, rows;
           int peermaster;
           while(numOfprocessor <= size){
                 //printf("I am here\n");
                 //============================================peer master==========================================//
                 if(rank%numOfprocessor==0){
                         //printf("I am rank %d\n", rank );

                         int avgrow = matrix_dimension_size / numOfprocessor;
                         int extraRow = matrix_dimension_size - numOfprocessor*avgrow;
                         offsetRow = (rank%numOfprocessor < extraRow) ? (avgrow+1) : avgrow;
                         copyMatrix(r[0],result,matrix_dimension_size);
                         int tmp;
                         tmp = offsetRow;
                         //=====================receive data from peer slave=======================================//
                         MPI_Recv(&r[1][0],matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE, rank+ previousNumOfProcess,rank+101,MPI_COMM_WORLD,&status);

                         // ====================send Data to peer slave==========================================//
                         int dest;
                         for( dest = 1; dest < numOfprocessor;dest++){
                                rows = (dest < extraRow) ? avgrow +1 : avgrow;
                                MPI_Send(&offsetRow, 1, MPI_INT, dest+rank, dest+rank+100, MPI_COMM_WORLD);
                                MPI_Send(&rows, 1, MPI_INT, dest+rank, dest+rank+100, MPI_COMM_WORLD);
                                MPI_Send(&r[0][offsetRow*matrix_dimension_size+ 0],rows*matrix_dimension_size,MPI_DOUBLE,rank+dest,rank+dest+100, MPI_COMM_WORLD);
                                MPI_Send(&r[1][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE,rank+dest,rank+dest+100,MPI_COMM_WORLD );
                                offsetRow += rows;
                         }

                         // ====================Calcluate on peer master========================================//
                         int k,i,j;
                         for (k=0; k<matrix_dimension_size; k++){
                                       for (i=0; i<tmp; i++) {
                                            result[i*matrix_dimension_size+k] = 0.0;
                                            for (j=0; j<matrix_dimension_size; j++)
                                                result[i*matrix_dimension_size+k] = result[i*matrix_dimension_size+k] + r[0][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
                                        }
                          }


                         //=====================Receive again from peer slave==================================//
                        int dest2;
                        for (dest2=1; dest2 < numOfprocessor; dest2++) {
                             int source = dest2+rank;
                             MPI_Recv(&offsetRow, 1, MPI_INT, source, rank+100,MPI_COMM_WORLD, &status); 
                             MPI_Recv(&rows, 1, MPI_INT, source, rank+100, MPI_COMM_WORLD, &status);
                             MPI_Recv(&result[offsetRow*matrix_dimension_size+0],rows*matrix_dimension_size, MPI_DOUBLE, source,rank+100,MPI_COMM_WORLD,&status);
                        }
                        //printf("I am rank %d\n", rank );

                 }
                 //==========================================peer slave=============================================//
                 else{
                         peermaster = rank;
                         while(peermaster % numOfprocessor !=0){
                            peermaster--;
                         }
                        // peermaster = (rank- previousNumOfProcess)/ numOfprocessor;
                         //===================== for some specific node, send to peer master =======================//
                         if(rank % previousNumOfProcess == 0 ){
                               //printf("I am rank %d\n", rank );
                               //printf("My master is  %d\n", peermaster);

                               MPI_Send(&result[0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE, peermaster, peermaster + 101, MPI_COMM_WORLD );
                         }

                         //=======================Receive from peer masrer==========================================//
                         
                         MPI_Recv(&offsetRow,1,MPI_INT,peermaster, rank+100, MPI_COMM_WORLD, &status);
                         MPI_Recv(&rows,1,MPI_INT,peermaster, rank+100, MPI_COMM_WORLD, &status);
                         MPI_Recv(&r[0][0], rows*matrix_dimension_size, MPI_DOUBLE,peermaster, rank  + 100, MPI_COMM_WORLD, &status);
                         MPI_Recv(&r[1][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE,peermaster, rank  + 100, MPI_COMM_WORLD, &status);

                         //=======================Calculate the result ============================================//
                         int k,i,j;
                         for (k=0; k<matrix_dimension_size; k++){
                                       for (i=0; i<rows; i++) {
                                            result[i*matrix_dimension_size+k] = 0.0;
                                            for (j=0; j<matrix_dimension_size; j++)
                                                result[i*matrix_dimension_size+k] = result[i*matrix_dimension_size+k] + r[0][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
                                        }
                          }

                         //=======================send back to peer  master ==========================================// 
                         MPI_Send(&offsetRow,1, MPI_INT, peermaster, peermaster+100,MPI_COMM_WORLD); 
                         MPI_Send(&rows, 1, MPI_INT, peermaster,peermaster+100,MPI_COMM_WORLD); 
                         MPI_Send(&result[0],rows*matrix_dimension_size,MPI_DOUBLE,peermaster,peermaster+100, MPI_COMM_WORLD);
                        // printf("I am rank %d\n", rank );

                 }
            MPI_Barrier(MPI_COMM_WORLD);
            previousNumOfProcess = numOfprocessor;
            numOfprocessor *= 2;
            //printf("number of processsor : %d \n", numOfprocessor);
           }

  }
  //=================================Number of processor is enough =====================================//
  else{

  }



   //=======================See result==================================================================//
  if(rank == 0){
      if (debug_perf == 0) {
        // print each of the sub matrices
        //for (i = 0; i < num_arg_matrices; ++i) {
          //  printf("argument matrix %d\n", i);
            //  print_matrix(r[i], matrix_dimension_size);
        //}
        printf("result matrix\n");
        print_matrix(result, matrix_dimension_size);
      } 
      /*else {
        double sum = 0.0;

        for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
            sum += result[resN][i];
        }
        printf("%f\n", sum);
      }
      */
  }
  

  MPI_Finalize();
}