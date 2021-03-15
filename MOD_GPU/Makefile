# Compiler 
CC = gcc
NVCC = nvcc
CFLAGS = -Wall -O2
NVCCFLAGS = -arch=sm_30 -O2 --compiler-options -Wall

# Include
INCLUDE += -I/home/x_fedlu/gsl/include

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))
CUSRC = $(wildcard *.cu)
CUOBJS = $(patsubst %.cu, %.o, $(CUSRC))

# Name of the executable
EXECUTABLE = stls

# Compile
all: $(EXECUTABLE)

%.o: %.c
	 $(CC) $(CFLAGS) $(INCLUDE) -c $<

%.o: %.cu
	 $(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $<

# Link
$(EXECUTABLE): $(OBJS) $(CUOBJS)
	 $(NVCC) $(LIB) $^ -o $@ -lgsl -lgslcblas -lm

clean:
	 rm *.o ${EXECUTABLE}
