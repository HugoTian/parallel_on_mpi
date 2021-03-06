#include <stdio.h>
#include <stdlib.h>
#include "gen_matrix.h"
#include "my_malloc.h"
#include "mpi.h"

int isPowerOfTwo(int n){
     while(n>1){
         if(n%2) return 0;
         n = n / 2;
     }
     return 1;
}
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
  int * info;

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

  info =(int *)my_malloc(sizeof(int) * size);
  
  if(num_arg_matrices == 1 && rank==0 ){
       result = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
        if (gen_sub_matrix(0, test_set, 0, result, 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
            printf("inconsistency in gen_sub_matrix\n");
            exit(1);
        }
        if (debug_perf == 0) {
        // print each of the sub matrices
    
            printf("argument matrix %d\n", i);
            print_matrix(result, matrix_dimension_size);
          
            printf("result matrix\n");
            print_matrix(result, matrix_dimension_size);
            //printf("The running time %f\n", end - start );
        } 
        else {
            double sum = 0.0;
            for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
            sum += result[i];
        }
            printf("%f\n", sum);
        }
        my_free(result);
        MPI_Finalize();
        return 0;
  }
  
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
           int fakesize = size;
           while(!isPowerOfTwo(fakesize)){
               fakesize++;
           }
           while(numOfprocessor <= fakesize){
                 //printf("I am here\n");
                 //============================================peer master==========================================//
                 if(rank%numOfprocessor==0){
                         //printf("I am rank %d\n", rank );
                         if(rank+numOfprocessor < size){
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
                          }else{
                                int newNumOfProcessor= size - rank;
                                int avgrow = matrix_dimension_size / newNumOfProcessor;
                                int extraRow = matrix_dimension_size - newNumOfProcessor*avgrow;
                                offsetRow = (rank% newNumOfProcessor < extraRow) ? (avgrow+1) : avgrow;
                                copyMatrix(r[0],result,matrix_dimension_size);
                                int tmp;
                                tmp = offsetRow;
                                //=====================receive data from peer slave=======================================//
                                if(rank + previousNumOfProcess < size){
                                    MPI_Recv(&r[1][0],matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE, rank+ previousNumOfProcess,rank+101,MPI_COMM_WORLD,&status);
                                }
                                // ====================send Data to peer slave==========================================//
                                int dest;
                                for( dest = 1; dest < newNumOfProcessor;dest++){
                                    rows = (dest < extraRow) ? avgrow +1 : avgrow;
                                    MPI_Send(&offsetRow, 1, MPI_INT, dest+rank, dest+rank+100, MPI_COMM_WORLD);
                                    MPI_Send(&rows, 1, MPI_INT, dest+rank, dest+rank+100, MPI_COMM_WORLD);
                                    MPI_Send(&r[0][offsetRow*matrix_dimension_size+ 0],rows*matrix_dimension_size,MPI_DOUBLE,rank+dest,rank+dest+100, MPI_COMM_WORLD);
                                    MPI_Send(&r[1][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE,rank+dest,rank+dest+100,MPI_COMM_WORLD );
                                    offsetRow += rows;
                                }

                              // ====================Calcluate on peer master========================================//
                              if(newNumOfProcessor !=1){
                                  int k,i,j;
                                  for (k=0; k<matrix_dimension_size; k++){
                                       for (i=0; i<tmp; i++) {
                                            result[i*matrix_dimension_size+k] = 0.0;
                                            for (j=0; j<matrix_dimension_size; j++)
                                                result[i*matrix_dimension_size+k] = result[i*matrix_dimension_size+k] + r[0][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
                                        }
                                  }
                              }
                              //=====================Receive again from peer slave==================================//
                              int dest2;
                              for (dest2=1; dest2 < newNumOfProcessor; dest2++) {
                                  int source = dest2+rank;
                                  MPI_Recv(&offsetRow, 1, MPI_INT, source, rank+100,MPI_COMM_WORLD, &status); 
                                  MPI_Recv(&rows, 1, MPI_INT, source, rank+100, MPI_COMM_WORLD, &status);
                                  MPI_Recv(&result[offsetRow*matrix_dimension_size+0],rows*matrix_dimension_size, MPI_DOUBLE, source,rank+100,MPI_COMM_WORLD,&status);
                              }
                          }
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
           //printf("fist reach here\n");
           int avgProcess = size / (num_arg_matrices/2);
           int extraProcess = size % (num_arg_matrices/2);
           int extraMatrix = num_arg_matrices % 2;
           int lastPeerProcess =rank;
           //double * num ;
           //num = (double *) my_malloc(sizeof(double) * size);

           //=============================Allocate Memory===============================================//
            r = (double **)my_malloc(sizeof(double *) * 2);
            r[0] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
            r[1] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
            result = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
           
            
            //===========================get the info of all process ===================================//
           
            int just;
            for(just = 0; just < size ; just++){
              info[just] = 0;
            }

           // ==============================Calculate Result on Each Process=============================================//
          int p = 0;
          int offsetNp = 0;
          int np ;
          for ( p=0; p< num_arg_matrices/2; p ++){
                np = ( p < extraProcess ) ? avgProcess+1 : avgProcess;
                info[offsetNp] =1;
                //=================use number of process with in range to calculate matrix p and matrix p and 1;
                if( rank < np+offsetNp && rank >= offsetNp){                     
                      int averow = matrix_dimension_size / np;
                      int extra = matrix_dimension_size % np;
                      int offset ,dest, rows, source;
                      //=============================allocate common memory=====================//
                      
                      if (gen_sub_matrix(0, test_set, 2*p+1, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                            printf("inconsistency in gen_sub_matrix\n");
                            exit(1);
                      }
                      // ==============master =====================//
                      if(rank == offsetNp){
                        if (gen_sub_matrix(0, test_set, 2*p, r[0], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                            printf("inconsistency in gen_sub_matrix\n");
                            exit(1);
                        }
                        if(np != 1){
                            offset = (0 < extra) ? averow +1 : averow;
                            int tmp = offset;
                            for(dest = 1; dest < np ; dest++){
                                //printf("I am rank %d, and I have %d group number\n", rank, np);

                                rows = (dest < extra) ? averow +1 : averow;
                                //printf("I am rank %dsending %d rows to task %d \n",rank , offset,dest);
                                MPI_Send(&offset, 1, MPI_INT, dest+rank, dest+rank+100, MPI_COMM_WORLD);
                                MPI_Send(&rows, 1, MPI_INT, dest+rank, dest+rank+100, MPI_COMM_WORLD);
                                MPI_Send(&r[0][offset*matrix_dimension_size + 0],rows*matrix_dimension_size, MPI_DOUBLE, dest+rank, rank+dest+100, MPI_COMM_WORLD);
                                //MPI_Send(&r[1][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE, dest,mtype, MPI_COMM_WORLD); 
                                offset = offset + rows;
                            } 
                            int k,i,j;
                            for (k=0; k<matrix_dimension_size; k++){
                                  for (i=0; i<tmp; i++) {
                                        result[i*matrix_dimension_size+k] = 0.0;
                                        for (j=0; j<matrix_dimension_size; j++)
                                             result[i*matrix_dimension_size+k] = result[i*matrix_dimension_size+k] + r[0][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
                                  }
                            }
                            /* wait for results from all slaves tasks */ 
                            int i2;
                            for (i2=1; i2< np; i2++) {
                                source = rank+i2;
                                MPI_Recv(&offset, 1, MPI_INT, source, rank+100,MPI_COMM_WORLD, &status); 
                                MPI_Recv(&rows, 1, MPI_INT, source, rank+100, MPI_COMM_WORLD, &status);
                                MPI_Recv(&result[offset*matrix_dimension_size+0],rows*matrix_dimension_size, MPI_DOUBLE, source,rank+100,MPI_COMM_WORLD,&status);
                                //printf("I am rank %d, and I am waiting  for %d\n", rank, offset);
                            }
                          }else{
                              mm(result,r[0],r[1],matrix_dimension_size);
                          }
                      }
                      else //=================slave ================//
                      {
                            //printf("I am rank %d, and my master is %d\n", rank,offsetNp );
                            MPI_Recv(&offset,1,MPI_DOUBLE,offsetNp, rank+100, MPI_COMM_WORLD, &status); 
                            MPI_Recv(&rows,1, MPI_INT, offsetNp, rank+100, MPI_COMM_WORLD, &status);
                            MPI_Recv(&r[0][0],rows*matrix_dimension_size,MPI_DOUBLE,offsetNp, rank+100,MPI_COMM_WORLD, &status);
                           // printf("I am rank %d, and I am waiting  for %d\n", rank, offsetNp);

                           
                            int k,i,j;
                            for (k=0; k<matrix_dimension_size; k++){
                                  for (i=0; i<rows; i++) {
                                        result[i*matrix_dimension_size+k] = 0.0;
                                        for (j=0; j<matrix_dimension_size; j++)
                                              result[i*matrix_dimension_size+k] = result[i*matrix_dimension_size+k] + r[0][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
                                  }
                            }

                             /*send back infomation to master*/
                             MPI_Send(&offset,1, MPI_INT, offsetNp, offsetNp+100,MPI_COMM_WORLD); 
                             MPI_Send(&rows, 1, MPI_INT, offsetNp, offsetNp+ 100,MPI_COMM_WORLD); 
                             MPI_Send(&result[0],rows*matrix_dimension_size,MPI_DOUBLE, offsetNp, offsetNp+100, MPI_COMM_WORLD);

                      }
                }
                offsetNp += np;
           }
            
        
           //printf("I am here\n");
           MPI_Barrier(MPI_COMM_WORLD);

           if(extraMatrix){
              if(rank == offsetNp - np){
                 copyMatrix(r[0], result,matrix_dimension_size);
                  if (gen_sub_matrix(0, test_set, num_arg_matrices-1, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                            printf("inconsistency in gen_sub_matrix\n");
                            exit(1);
                        }
                 mm(result,r[0],r[1],matrix_dimension_size);
              }
           }
           

           //===============================Combine result from bottom up ========================================//
           int just2 ,sumTmp;
           int previousNumOfProcess;
           sumTmp = 0;
           for(just2 = 0 ; just2 < size ; just2++){
                if(info[just2]){
                    if(sumTmp == 0)
                      sumTmp = 1;
                    else{
                      sumTmp = 0;
                      info[just2] = -1;
                    }
                }
           }
           int tmpsize = size;
           int numOfprocessor = 2;
           int offsetRow, rows;
           int peermaster;
           int sum ;
           sum= 0 ;
           int t11;
           for (t11 = 0 ; t11 < size ; t11++){
                if(info[t11])
                    sum++;
           }
          // printf("sum is %d\n", sum);
           MPI_Barrier(MPI_COMM_WORLD);
           //printf("I am here\n");

           
           while( sum!= 1){
                 //printf("I am here\n");
                  //if(rank == 0){
                  //    int tmptmp =0;
                  ////    for(tmptmp = 0 ; tmptmp < size ; tmptmp++){
                   //       printf("%d ", info[tmptmp] );
                   //   }
                 // }
                 // printf("\n");
                 //============================================peer master==========================================//
                 if(info[rank]==1){
                         //printf("I am rank %d\n", rank );
                         numOfprocessor = 1;
                         int tmpRank = rank+1 ;
                         while(info[tmpRank] != 1 && tmpRank < size){
                            numOfprocessor++;
                            tmpRank++;
                         }
                         previousNumOfProcess = 1;
                         int tmpRank2 = rank + 1;
                         while( tmpRank2 < size && info[tmpRank2] != -1 ){
                             previousNumOfProcess ++;
                             tmpRank2++;
                         }

                         if(tmpRank2 != size){
                              //printf("I am rank %d, and I have member %d", rank  , numOfprocessor);

                              int avgrow = matrix_dimension_size / numOfprocessor;
                              int extraRow = matrix_dimension_size  % numOfprocessor;
                              offsetRow = (0 < extraRow) ? (avgrow+1) : avgrow;
                              copyMatrix(r[0],result,matrix_dimension_size);
                              int tmp;
                              tmp = offsetRow;
                              //=====================receive data from peer slave=======================================//
                              //printf("I am rank %d, and I am waiting for %d \n", rank, rank+previousNumOfProcess );

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
                        
                 }
                 //==========================================peer slave=============================================//
                 else{
                         peermaster = rank;
                         while(info[peermaster] != 1){
                            peermaster--;
                         }

                         int working = peermaster+1;
                         while(working< size &&  info[working] != -1  ){
                              working++;
                         }
                         if(working != size){

                         lastPeerProcess = peermaster;
                        // printf("I am rank %d and my peer master is %d \n", rank, peermaster );

                        // peermaster = (rank- previousNumOfProcess)/ numOfprocessor;
                         //===================== for some specific node, send to peer master =======================//
                         if(info[rank] == -1 ){
                               //printf("I am rank %d\n", rank );
                               //printf("My master is  %d\n", peermaster);
                               //printf("I am rank %d and I send to peer master is %d \n", rank, peermaster );

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

                 }
            MPI_Barrier(MPI_COMM_WORLD);
            // update the info array
            //previousNumOfProcess = numOfprocessor;
            //numOfprocessor *= 2;
            //printf("number of processsor : %d \n", numOfprocessor);
            int t1 =0;
            int sum3 =0;
            for (t1 = 0 ; t1 < size ; t1++){
                if(info[t1] == 1){
                     if(sum3 == 0 )
                        sum3 = 1;
                     else{
                         sum3 = 0 ;
                         info[t1] = -1;
                     }
                }
                else if(info[t1] == -1){
                      info[t1] =0;
                }
            }

            sum= 0 ;
            for (t1 = 0 ; t1 < size ; t1++){
                if(info[t1] != 0 )
                    sum++;
            }
            MPI_Barrier(MPI_COMM_WORLD);


          }


  }


  end = MPI_Wtime();


   //=======================See result==================================================================//
  if(rank == 0){
      if (debug_perf == 0) {
        // print each of the sub matrices
        //for (i = 0; i < num_arg_matrices; ++i) {
        //    printf("argument matrix %d\n", i);
        //      print_matrix(r[i], matrix_dimension_size);
        //}
        printf("result matrix\n");
        
        print_matrix(result, matrix_dimension_size);


        printf("The running time %f\n", end - start );
      } 
      else {
        double sum = 0.0;

        for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
            sum += result[i];
        }
        printf("%f\n", sum);
      }
      
  }
  my_free(result);
  my_free(r[0]);
  my_free(r[1]);
  my_free(info);
  MPI_Finalize();

}