# Default type of integration scheme
ifndef $(ITG)
ITG = cquad
endif

# Name of output file
ifndef $(EXE)
EXE = stls
endif

# Path for the GSL library
ifndef $(GSL)
GSL = -I/home/x_fedlu/gsl/include
endif

# Run appropriate make files
.PHONY : all
all:
	cd $(ITG) && make $@ EXE=$(EXE) GSL=$(GSL) && mv $(EXE) ..


.PHONY : clean
clean:
	cd $(ITG) && make $@
	rm $(EXE)

