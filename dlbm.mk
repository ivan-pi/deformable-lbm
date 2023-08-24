

FC = gfortran
FFLAGS = -Wall -fcheck=all -fimplicit-none -Og -g

LDFLAGS = 

dlbm: dlbm.f90 dlbm_cmd.f90 dlbm_io.f90 deformable_lbm.f90 
	$(FC) $(FFLAGS) -o $@ $^

.phony: clean cleanobj cleanmod cleantxt

clean: cleanobj cleanmod cleantxt
	rm -f dlbm

cleanobj :
	rm -f *.o

cleanmod :
	rm -f *.mod *__genmod.f90

cleantxt :
	rm -f *.txt