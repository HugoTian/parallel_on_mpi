test_mm: mpi.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h
	mpicc -O3 -o mpi  mpi.c gen_matrix.c my_malloc.c

cilk_mm: cilk.c gen_matrix.c my_malloc.c  gen_matrix.h my_malloc.h
	icc -O3 -o cilk_mm cilk.c  gen_matrix.c my_malloc.c
run_debug:
	./test_mm 0 2 10

run_performance:
	./test_mm 1 0 10

clean:
	rm *~; rm *.exe
