OS=$(shell uname)
CC?=gcc
CFLAGS=-O3 -Wall -Wextra -std=c99 -march=native -ffast-math

# Numpy expects differents extensions for a dynamic library on Linux vs. OSX
# Also, the C compiler option to build a dynamic library is different
SONAME=libffa.so
SOOPT=-soname
ifeq ($(OS),Darwin)
	SONAME=libffa.dylib
	SOOPT=-install_name
endif

$(SONAME): libffa.o
	$(CC) $(CFLAGS) -shared -Wl,$(SOOPT),$(SONAME) -o $(SONAME) libffa.o

libffa.o: libffa.c *.h
	$(CC) $(CFLAGS) -c -fPIC libffa.c -o libffa.o

all: $(SONAME)

.PHONY: clean 

clean:
	rm -vf *.o
	rm -vf $(SONAME)
