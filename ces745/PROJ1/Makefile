default: tsp_serial tsp_openmp

tsp_serial: tsp_serial.c
	icc -O2 tsp_serial.c -o tsp_serial

tsp_openmp: tsp_openmp.c
	icc -O2 -openmp tsp_openmp.c -o tsp_openmp
