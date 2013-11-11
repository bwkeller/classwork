#!/bin/bash
while true
do
	min=$RANDOM
	max=$RANDOM
	let "max *=100"
	make CFLAGS="-DKMIN=$min -DKMAX=$max"
	parallel=`./parallelprimes | grep maxval | cut -d' ' -f 2`
	serial=`./primes | grep maxval | cut -d' ' -f 2`
	if [ $parallel -ne $serial ]
	then
		echo "FAILED parallel: "$parallel" serial: "$serial
		exit 1
	fi
	echo "PASSED parallel: "$parallel" serial: "$serial
	make clean
done
