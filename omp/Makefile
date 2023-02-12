######################################################################
# @author      : Evan Gerritz (evan.gerritz@yale.edu)
# @file        : Makefile
# @created     : Friday Dec 16, 2022 15:26:36 EST
######################################################################

IDIR =./include
CC=icc
CFLAGS=-I$(IDIR) -O3 -Wall -xHost -fno-alias -std=c99 -qopenmp

ODIR=obj

LIBS=

_OBJ = ga-omp-sort.o timing.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(OBJ)
	$(CC) -o ga-omp-sort $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o

