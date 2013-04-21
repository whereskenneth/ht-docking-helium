CC = g++
CFLAGS = -c
LDFLAGS = -L${SHARED_ROOT}/usr/lib -lopenbabel -lboost_program_options -I${SHARED_ROOT}/usr/include
GFLAGS=

all: sdftopdbqt

sdftopdbqt:  sdftopdbqt.o 
	$(CC) $(LDFLAGS)  sdftopdbqt.o -o sdftopdbqt

sdftopdbqt.o: sdftopdbqt.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) sdftopdbqt.cpp






	



clean:
	rm -rf sdftopdbqt.o sdftopdbqt
