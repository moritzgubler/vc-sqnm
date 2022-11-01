FC := gfortran

libs:= -llapack -lblas

ftn_flags := -g

program := sqnm_test.exe


ftn_src := $(wildcard *.f90)
ftn_obj := $(patsubst %.f90, %.o, $(ftn_src))


.PHONY:	clean

$(program): $(ftn_obj)
	$(FC) $(ftn_flags) $^ $(libs) -o $@

%.o: %.f90
	$(FC) $(ftn_flags) -c $< -o $@


clean:
	rm -rf *.o *.mod $(program)


#dependencies
sqnm.o: historylist.o
periodic_optimizer.o: sqnm.o
main.o: periodic_optimizer.o