# Compiler 
CC = gcc
CFLAGS = -Wall

# Include
INCLUDE += -I/home/x_fedlu/gsl/include

# Libraries
LIB = -static -L/home/x_fedlu/gsl/lib

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

# Name of the executable
EXECUTABLE = stls

# Compile
all: $(EXECUTABLE)

%.o: %.c
	 $(CC) $(CFLAGS) -fopenmp $(INCLUDE) -c $<

# Link
$(EXECUTABLE): $(OBJS)
	 $(CC) -fopenmp $(LIB) $^ -o $@ -lgsl -lgslcblas -lm

clean:
	 rm *.o ${EXECUTABLE}