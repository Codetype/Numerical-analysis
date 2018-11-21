loadpath:
	LD_LIBRARY_PATH=/root/gsl/lib
	export LD_LIBRARY_PATH
	
runzad1:
	gcc -o zad1 zad1.c -lm
	./zad1 0 10 100
	./zad1 0 10 1000
	./zad1 0 10 10000
	rm -f zad1

runzad2:
	gcc -Wall -I/root/gsl/include -c zad2.c
	gcc -L/root/gsl/lib zad2.o -lgsl -lgslcblas -lm
	./a.out 10 6 1 9 100
	./a.out 10 6 1 9 1000
	./a.out 10 6 1 9 10000
	rm -f zad2.o a.out