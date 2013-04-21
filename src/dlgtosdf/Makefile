CC = g++
CFLAGS = -c
LDFLAGS = -L${SHARED_ROOT}/usr/lib -lopenbabel -lboost_program_options -I${SHARED_ROOT}/usr/include
GFLAGS=

all: dlgtosdf

dlgtosdf:  dlgtosdf.o WriteToSDF.o
	$(CC) $(LDFLAGS)  WriteToSDF.o dlgtosdf.o -o dlgtosdf

dlgtosdf.o: dlgtosdf.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) dlgtosdf.cpp

WriteToSDF.o: WriteToSDF.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) WriteToSDF.cpp





	



clean:
	rm -rf WriteToSDF.o dlgtosdf.o dlgtosdf
