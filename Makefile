# Compiler 
CC = gcc
CFLAGS = -Wall -O2 -std=gnu99

# Include
INCLUDE += -I/home/flc/nobackup/gsl/include
LIB += -L/home/flc/nobackup/gsl/lib

# Name of output file
EXE = stls

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

# Compile
all: $(EXE)

%.o: %.c
	 $(CC) $(CFLAGS) -fopenmp $(INCLUDE) -c $<

# Link
$(EXE): $(OBJS)
	 $(CC) $(LIB) -fopenmp $^ -o $@ -lgsl -lgslcblas -lm

clean:
	 rm *.o $(EXE)
