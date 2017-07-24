CC=gcc
CFLAGS=-O3 -Wall -std=c99 --fast-math
VECFLAGS=-fopt-info -ftree-vectorizer-verbose=1
OMPFLAGS=-fopenmp

libffa.so: libffa.o
	$(CC) $(CFLAGS) $(OMPFLAGS) -shared -Wl,-soname,libffa.so -o libffa.so libffa.o

libffa.o: libffa.c
	$(CC) $(CFLAGS) $(OMPFLAGS) $(VECFLAGS) -c -fPIC libffa.c -o libffa.o

all: libffa.so

.PHONY: clean

clean:
	rm -vf *.o
	rm -vf *.so
