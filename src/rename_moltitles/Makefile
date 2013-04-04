CC = g++
CFLAGS = -c
LDFLAGS = -L/Users/ktussey/usr/local/lib -lopenbabel -lboost_program_options -I/Users/ktussey/usr/local/include
GFLAGS=

all: rename_moltitles

rename_moltitles:  rename_moltitles.o 
	$(CC) $(LDFLAGS)  rename_moltitles.o -o rename_moltitles

rename_moltitles.o: rename_moltitles.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) rename_moltitles.cpp






	



clean:
	rm -rf rename_moltitles.o rename_moltitles
