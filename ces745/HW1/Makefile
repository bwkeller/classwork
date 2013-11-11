default: parallelprimes primes

parallelprimes: parallelprimes.c
	cc $(CFLAGS) parallelprimes.c -lm -openmp -o parallelprimes

primes: primes.c
	cc $(CFLAGS) primes.c -lm -o primes
clean:
	@rm parallelprimes primes
