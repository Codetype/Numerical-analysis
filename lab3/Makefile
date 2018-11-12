run: runzad1 clean

loadpath:
	LD_LIBRARY_PATH=/root/gsl/lib
	export LD_LIBRARY_PATH

runzad1:
	gcc -Wall -I/root/gsl/include -c zad1.c
	gcc -L/root/gsl/lib zad1.o -lgsl -lgslcblas -lm
	./a.out

runzad2:
	gcc -o zad2 zad2.c
	./zad2

runzad3:
	gcc -o zad3 zad3.c
	./zad3

runzad4:
	gcc -Wall -I/root/gsl/include -c zad4.c
	gcc -L/root/gsl/lib zad4.o -lgsl -lgslcblas -lm
	./a.out

clean:
	rm -f *.o *.out zad2 zad3
