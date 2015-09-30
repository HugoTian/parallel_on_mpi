make test_mm
make mpi_mm
echo "run test set $1 with size $2 on $3 number of process"
./test_mm 0 $1 $2 > result1.txt
mpirun -np $3 ./mult_mpi 0 $1 $2 > result2.txt
diff result1.txt result2.txt
