CC = g++
CFLAGS = -c
LDFLAGS = -L/Users/ktussey/usr/local/lib -lopenbabel -lboost_program_options -I/Users/ktussey/usr/local/include
GFLAGS=

all: sdftopdbqt

sdftopdbqt:  sdftopdbqt.o 
	$(CC) $(LDFLAGS)  sdftopdbqt.o -o sdftopdbqt

sdftopdbqt.o: sdftopdbqt.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) sdftopdbqt.cpp






	



clean:
	rm -rf sdftopdbqt.o sdftopdbqt
