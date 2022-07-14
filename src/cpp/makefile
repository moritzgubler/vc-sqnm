CPPC := g++
FC := gfortran

libs:= -llapack -lblas -lgfortran

cpp_flags := -g -Wall
ftn_flags := -g -Wall

program := sqnm_test.exe


cpp_src := $(wildcard *.cpp)
ftn_src := $(wildcard *.f90)

cpp_obj := $(patsubst %.cpp, %.o, $(cpp_src))
ftn_obj := $(patsubst %.f90, %.ftn_o, $(ftn_src))


.PHONY:	clean

$(program): $(cpp_obj) $(ftn_obj)
	$(CPPC) $(cpp_flags) $^ $(libs) -o $@

%.ftn_o: %.f90
	$(FC) $(ftn_flags) -c $< -o $@

%.o: %.cpp
	$(CPPC) -c $(cpp_flags) $< -o $@


clean:
	rm -rf *.o *.ftn_o $(program)