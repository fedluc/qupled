# Default type of integration scheme
ifndef $(ITG)
ITG = cquad
endif

# Name of output file
ifndef $(EXE)
EXE = stls
endif

# Run appropriate make files
.PHONY : all
all:
	cd $(ITG) && make $@ EXE=$(EXE) && mv $(EXE) ..


.PHONY : clean
clean:
	cd $(ITG) && make $@
	rm $(EXE)

