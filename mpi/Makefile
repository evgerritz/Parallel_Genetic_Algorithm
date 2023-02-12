#
# This Makefile assumes that the following module files are loaded:
#
#   iomkl
#

# Set the compiler command
CC = mpicc

# Set the compiler options
CFLAGS = -g3 -O3 -std=c99 # -xHost -fno-alias

EXECS = ga-island
OBJS = ga-island.o

all: $(EXECS)

ga-island: ga-island.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

# Be careful not to remove rwork.o file
clean:
	rm -f $(EXECS) $(OBJS)
