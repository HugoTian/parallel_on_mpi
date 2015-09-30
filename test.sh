make test_mm
make sec_mpi
echo "Run test set $1 and matrix size $2 in $3 processes"
./test_mm 0 $1 $2 > res1.txt
mpirun -np $3 ./mpi2 0 $1 $2 > res2.txt
diff res1.txt res2.txt
