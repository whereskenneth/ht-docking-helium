CC = g++
CFLAGS = -c
LDFLAGS = -L${SHARED_ROOT}/usr/lib -lboost_program_options -I${SHARED_ROOT}/usr/include
GFLAGS=

all: ad4statistics

ad4statistics:  ad4statistics.o statistics_functions.o write.o
	$(CC)	$(LDFLAGS)	ad4statistics.o statistics_functions.o write.o -o ad4statistics

ad4statistics.o:  ad4statistics.cpp
	$(CC)	$(CFLAGS)	$(LDFLAGS)	ad4statistics.cpp

statistics_functions.o: statistics_functions.cpp
	$(CC)	$(CFLAGS)	$(LDFLAGS) 	statistics_functions.cpp

write.o:	write.cpp
	$(CC)	$(CFLAGS)	$(LDFLAGS)	write.cpp





	



clean:
	rm ad4statistics ad4statistics.o statistics_functions.o write.o
