CC=gcc
CFLAGS=-O3 -Wall -std=c99 --fast-math

libffa.so: libffa.o
	$(CC) $(CFLAGS) -shared -Wl,-soname,libffa.so -o libffa.so libffa.o

libffa.o: libffa.c
	$(CC) $(CFLAGS) -c -fPIC libffa.c -o libffa.o

all: libffa.so

.PHONY: clean

clean:
	rm -vf *.o
