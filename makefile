OBJS = clixo.o
CC = g++
CFLAGS11 = -Wall -O4 -c -std=c++11 -fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations -I .
LFLAGS11 = -Wall -O4 -std=c++11 -fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations -I .

all: clixo 

clixo: $(OBJS)
	$(CC) $(LFLAGS11) clixo.o -o clixo

clixo.o: clixo.cpp util.h dag.h graph_undirected.h dagConstruct.h graph_undirected_bitset.h nodeDistanceObject.h
	$(CC) $(CFLAGS11) clixo.cpp

clean:
	rm *.o
