make test_mm
make sec_mpi

for i in `seq 1 10`; do
	for j in `seq 1 10`; do
		./test_mm 0 3 $i > res1.txt
		mpirun -np $j ./mpi2 0 3 $i > res2.txt
		diff res1.txt res2.txt
		if [ $? -eq 0 ]  ; then
			echo "Run test set 3 and matrix size $i in $j processes success"
		else
			echo "Fail for size $i and process $j"
		fi
		
	done
done
