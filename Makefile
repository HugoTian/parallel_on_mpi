test_mm: test_mm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h
	gcc -g -DDEBUG test_mm.c gen_matrix.c my_malloc.c -o test_mm
my_mm:  matrix_multi.cilk  gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h
	cilkc -D_XOPEN_SOURCE=600 -D_POSIX_C_SOURCE=200809L  gen_matrix.c my_malloc.c matrix_multi.cilk -o mult

mpi_mm: matrix_mpi.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h
	mpicc matrix_mpi.c gen_matrix.c my_malloc.c -o mult_mpi
run_debug:
	./test_mm 0 0 10

run_performance:
	./test_mm 1 0 10

run_my:
	./mult 0 0 10
clean:
	rm *~; rm *.exe
