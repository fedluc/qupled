# Compiler 
CC = gcc
CFLAGS = -Wall -O2

# Include
INCLUDE += -I/home/x_fedlu/gsl/include

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

# Name of the executable
EXECUTABLE = stls

# Compile
all: $(EXECUTABLE)

%.o: %.c
	 $(CC) $(CFLAGS) $(INCLUDE) -c $<

# Link
$(EXECUTABLE): $(OBJS)
	 $(CC) $^ -o $@ -lgsl -lgslcblas -lm

clean:
	 rm *.o ${EXECUTABLE}
