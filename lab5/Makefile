loadpath:
	LD_LIBRARY_PATH=/root/gsl/lib
	export LD_LIBRARY_PATH

runzad1:
	gcc -Wall -I/root/gsl/include -c zad1.c
	gcc -L/root/gsl/lib zad1.o -lgsl -lgslcblas -lm
	./a.out 4 0.00001
	rm -f zad1.o a.out

runzad2:
	gcc -Wall -I/root/gsl/include -c zad2.c
	gcc -L/root/gsl/lib zad2.o -lgsl -lgslcblas -lm
	./a.out 4 0.00001
	rm -f zad2.o a.out

runzad3:
	gcc -Wall -I/root/gsl/include -c zad3.c
	gcc -L/root/gsl/lib zad3.o -lgsl -lgslcblas -lm
	./a.out 4 0.00001
	rm -f zad3.o a.out

runzad4:
	gcc -Wall -I/root/gsl/include -c zad4.c
	gcc -L/root/gsl/lib zad4.o -lgsl -lgslcblas -lm
	./a.out 10 10
	rm -f zad4.o a.out

run_testiternum_zad4:
	gcc -Wall -I/root/gsl/include -c zad4.c
	gcc -L/root/gsl/lib zad4.o -lgsl -lgslcblas -lm
	./a.out 5 10
	./a.out 5 20
	./a.out 5 30
	./a.out 5 40
	./a.out 5 50
	./a.out 5 60
	./a.out 5 70
	./a.out 5 1000
	./a.out 5 10000
	rm -f zad4.o a.out


run_testsize_zad4:
	gcc -Wall -I/root/gsl/include -c zad4.c
	gcc -L/root/gsl/lib zad4.o -lgsl -lgslcblas -lm
	./a.out 5 100
	./a.out 10 100
	./a.out 15 100
	./a.out 20 100
	./a.out 25 100
	rm -f zad4.o a.out