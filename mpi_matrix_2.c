#include <stdio.h>
#include <stdlib.h>
#include "gen_matrix.h"
#include "my_malloc.h"
#include "mpi.h"

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


  int averow = matrix_dimension_size / size;
  int extra = matrix_dimension_size % size;
  
  int offset, mtype ,dest, rows, source;
  int from_master = 111;
  int from_slave = 1111;
  int Master = 0;

  //====================== allocate memory==========================================================//
  if(rank == Master){
  	 r = (double **)my_malloc(sizeof(double *) * 2);
  	 //for (i = 0; i < num_arg_matrices; ++i) {
    	 r[0] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
    	 r[1] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
    	 if (gen_sub_matrix(0, test_set, 0, r[0], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
     		 printf("inconsistency in gen_sub_matrix\n");
     		 exit(1);
    	 }
  	 //}  
  	 result = (double **)my_malloc(sizeof(double *) * 2);
     result[0] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
     result[1] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
  }
  else{
  	rp = (double **)my_malloc(sizeof(double *) * 3);
  	rp[0] = (double *)my_malloc(sizeof(double) * (averow+1) * matrix_dimension_size);
  	rp[1] = (double *)my_malloc(sizeof(double) * matrix_dimension_size * matrix_dimension_size);
  	rp[2] = (double *)my_malloc(sizeof(double) * (averow+1) * matrix_dimension_size);
  }


  /* 
  *  Special case calcute result 0 first
  *  Matrix Multiplication
  *  result[0] = r[0] X r[1]
  */

  //====================Master task ======================================================================//

  if(rank == Master){
  		//printf("into master job\n");
  		dest = 0;
  		offset = (dest < extra) ? averow +1 : averow;
  		int tmp = offset;
  		mtype = from_master;
  		if(size !=1){
  			//=====send task to slave ==========================================//
  			for(dest = 1; dest < size ; dest++){
  				rows = (dest < extra) ? averow +1 : averow;
  				//printf("sending %d rows to task %d \n",rows,dest);
				MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&r[0][offset*matrix_dimension_size + 0],rows*matrix_dimension_size, MPI_DOUBLE, dest,mtype, MPI_COMM_WORLD);
				//MPI_Send(&r[1][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE, dest,mtype, MPI_COMM_WORLD); 
				offset = offset + rows;
  			}	
  			//=============================task for master=====================//
  			if (gen_sub_matrix(0, test_set, 1, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
     		   printf("inconsistency in gen_sub_matrix\n");
     		   exit(1);
    	    }
  			int k,i,j;
	    	for (k=0; k<matrix_dimension_size; k++){

				for (i=0; i<tmp; i++) {
					 result[0][i*matrix_dimension_size+k] = 0.0;
					 for (j=0; j<matrix_dimension_size; j++)
						  result[0][i*matrix_dimension_size+k] = result[0][i*matrix_dimension_size+k] + r[0][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
				}
			}
  			/* wait for results from all slaves tasks */ 
  			mtype = from_slave;
  			for (i=1; i< size; i++) {
				source = i;
				MPI_Recv(&offset, 1, MPI_INT, source, mtype,MPI_COMM_WORLD, &status); 
				MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
				MPI_Recv(&result[0][offset*matrix_dimension_size+0],rows*matrix_dimension_size, MPI_DOUBLE, source,mtype,MPI_COMM_WORLD,&status);
			}
		}
		else{
      if (gen_sub_matrix(0, test_set, 1, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
           printf("inconsistency in gen_sub_matrix\n");
           exit(1);
      }
			mm(result[0],r[0],r[1],matrix_dimension_size);
    }
    
  }

 //====================Slave task=========================================================================//
  if (rank > Master) {
		mtype = from_master;
		MPI_Recv(&offset,1,MPI_DOUBLE,Master, mtype, MPI_COMM_WORLD, &status); 
		MPI_Recv(&rows,1, MPI_INT, Master, mtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&rp[0][0],rows*matrix_dimension_size,MPI_DOUBLE,Master,mtype,MPI_COMM_WORLD, &status);
	   // MPI_Recv(&rp[1][0],matrix_dimension_size*matrix_dimension_size,MPI_DOUBLE, Master,mtype,MPI_COMM_WORLD,&status);
        
        //if (gen_sub_matrix(0, test_set, 0, rp[0], offset, offset + rows -1 , 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
     	//	 printf("inconsistency in gen_sub_matrix\n");
     	//	 exit(1);
    	// }
	    if (gen_sub_matrix(0, test_set, 1, rp[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
     		 printf("inconsistency in gen_sub_matrix\n");
     		 exit(1);
    	 }
	    int k,i,j;
	    for (k=0; k<matrix_dimension_size; k++){
			for (i=0; i<rows; i++) {
				 rp[2][i*matrix_dimension_size+k] = 0.0;
				for (j=0; j<matrix_dimension_size; j++)
					rp[2][i*matrix_dimension_size+k] = rp[2][i*matrix_dimension_size+k] + rp[0][i*matrix_dimension_size+j] * rp[1][j*matrix_dimension_size+k]; 
			}
		}

		/*send back infomation to master*/
		mtype = from_slave;
		MPI_Send(&offset,1, MPI_INT, Master, mtype,MPI_COMM_WORLD); 
		MPI_Send(&rows, 1, MPI_INT, Master,mtype,MPI_COMM_WORLD); 
		MPI_Send(&rp[2][0],rows*matrix_dimension_size,MPI_DOUBLE,Master,mtype, MPI_COMM_WORLD);

  }
  MPI_Barrier(MPI_COMM_WORLD);
  //====================USE loop to calculate result in cascade =======================================//
  //if(rank == 0){
  //printf("result matrix\n");
  //print_matrix(result[0], matrix_dimension_size);
  //}
  int number = 2;
  int resN = 0; 
  
  for (number = 2; number < num_arg_matrices ; number++ ){
  		//offset = 0;

  	    //====================Master task ======================================================================//

  		if(rank == Master){
  				//printf("into master job\n");
  				dest = 0;
  				offset = (dest < extra) ? averow +1 : averow;
  				int tmp = offset;
  				mtype = from_master;
  				if(size !=1){
  					//=====send task to slave ==========================================//
  					for(dest = 1; dest < size ; dest++){
  						rows = (dest < extra) ? averow +1 : averow;
  						//printf("sending %d rows to task %d \n",rows,dest);
						MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
						MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
						MPI_Send(&result[resN][offset*matrix_dimension_size + 0],rows*matrix_dimension_size, MPI_DOUBLE, dest,mtype, MPI_COMM_WORLD);
						//MPI_Send(&r[number][0], matrix_dimension_size*matrix_dimension_size, MPI_DOUBLE, dest,mtype, MPI_COMM_WORLD); 
						offset = offset + rows;
  					}	
  					//=============================task for master=====================//
  					if (gen_sub_matrix(0, test_set, number, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
     		   				printf("inconsistency in gen_sub_matrix\n");
     		   				exit(1);
    	    		}
  					int k,i,j;
	    			for (k=0; k<matrix_dimension_size; k++){
						for (i=0; i<tmp; i++) {
					 		result[resN ^ 0x1][i*matrix_dimension_size+k] = 0.0;
					 		for (j=0; j<matrix_dimension_size; j++)
						  		result[resN ^ 0x1][i*matrix_dimension_size+k] = result[resN ^ 0x1][i*matrix_dimension_size+k] + result[resN][i*matrix_dimension_size+j] * r[1][j*matrix_dimension_size+k]; 
						}
					}
  					// wait for results from all slaves tasks 
  					mtype = from_slave;
  					for (i=1; i< size; i++) {
						source = i;
						MPI_Recv(&offset, 1, MPI_INT, source, mtype,MPI_COMM_WORLD, &status); 
						MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
						MPI_Recv(&result[resN ^ 0x1][offset*matrix_dimension_size+0],rows*matrix_dimension_size, MPI_DOUBLE, source,mtype,MPI_COMM_WORLD,&status);
					}
				}
				else{
          if (gen_sub_matrix(0, test_set, number, r[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
                  printf("inconsistency in gen_sub_matrix\n");
                  exit(1);
          }
					mm(result[resN ^ 0x1],result[resN],r[1],matrix_dimension_size);
        }
  	    }

 		//====================Slave task=========================================================================//
  		if (rank > Master) {
				mtype = from_master;
				MPI_Recv(&offset,1,MPI_DOUBLE,Master, mtype, MPI_COMM_WORLD, &status); 
				MPI_Recv(&rows,1, MPI_INT, Master, mtype, MPI_COMM_WORLD, &status);
	    		MPI_Recv(&rp[0][0],rows*matrix_dimension_size,MPI_DOUBLE,Master,mtype,MPI_COMM_WORLD, &status);
	   			//MPI_Recv(&rp[1][0],matrix_dimension_size*matrix_dimension_size,MPI_DOUBLE, Master,mtype,MPI_COMM_WORLD,&status);
	   			if (gen_sub_matrix(0, test_set, number, rp[1], 0, matrix_dimension_size - 1, 1, 0, matrix_dimension_size - 1, 1, 1) == NULL) {
     		 		printf("inconsistency in gen_sub_matrix\n");
     		 		exit(1);
    	 		}
	   			int k,i,j;
	    		for (k=0; k<matrix_dimension_size; k++){
					for (i=0; i<rows; i++) {
				 		rp[2][i*matrix_dimension_size+k] = 0.0;
						for (j=0; j<matrix_dimension_size; j++)
							rp[2][i*matrix_dimension_size+k] = rp[2][i*matrix_dimension_size+k] + rp[0][i*matrix_dimension_size+j] * rp[1][j*matrix_dimension_size+k]; 
					}
				}

				//send back infomation to master
				mtype = from_slave;
				MPI_Send(&offset,1, MPI_INT, Master, mtype,MPI_COMM_WORLD); 
				MPI_Send(&rows, 1, MPI_INT, Master,mtype,MPI_COMM_WORLD); 
				MPI_Send(&rp[2][0],rows*matrix_dimension_size,MPI_DOUBLE,Master,mtype, MPI_COMM_WORLD);

 		}

 		MPI_Barrier(MPI_COMM_WORLD);
 		resN = resN ^ 0x1;
  }


  end = MPI_Wtime();
  MPI_Finalize();

  //=======================See result==================================================================//
  if(rank == 0){
   		if (debug_perf == 0) {
    		// print each of the sub matrices
    		//for (i = 0; i < num_arg_matrices; ++i) {
     	    //	printf("argument matrix %d\n", i);
      	    //	print_matrix(r[i], matrix_dimension_size);
    	 	//}
    	 	printf("result matrix\n");
    	 	print_matrix(result[resN], matrix_dimension_size);
        printf("The program Running time for mpi2 is %f\n", end -start );
   		} else {
    		double sum = 0.0;

    		for (i = 0; i < matrix_dimension_size * matrix_dimension_size; ++i) {
      			sum += result[resN][i];
    		}
    		printf("%f\n", sum);
  		}
  }
  

}
