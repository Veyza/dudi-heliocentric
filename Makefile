# This file is a part of DUDI-heliocentric, the Fortran-95 implementation 
# of the two-body model for the dynamics of dust ejected from an atmosphereless
# body moving around the Sun
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

objs = const.o define_types.o help.o distributions_fun.o data_in.o data_out.o twobody_fun.o DUDIhc.o
f90s = const.f90 define_types.f90 help.f90 distributions_fun.f90 data_in.f90 data_out.f90 gu.f90 twobody_fun.f90 DUDIhc.f90

example_image : example
	./dudihc
	python3 show_image.py

phaethon : phaethon_compile
	./phaethon_dudihc
	python3 plot_Fig10.py

example : example.o $(objs)
	gfortran -fopenmp $(objs) example.o -o dudihc

phaethon_compile : $(objs) phaethon_input.o phaethon.o
	gfortran -fopenmp -o phaethon_dudihc $(objs) phaethon_input.o phaethon.o

select_method : select_method_comp ./input_data_files/orbit_and_time_test.dat
	./select_method_dudihc

select_method_comp : select_method.o  $(objs)
	gfortran -fopenmp -o select_method_dudihc $^

select_method.o : select_method.f90 $(objs)
	gfortran -c $< -fopenmp

example.o : example.f90 $(objs)
	gfortran -c $< -fopenmp

phaethon.o : phaethon.f90 $(objs) phaethon_input.o
	gfortran -c $< -fopenmp
	
phaethon_input.o : phaethon_input.f90 help.o const.o
	gfortran -c $< -fopenmp

DUDIhc.o : DUDIhc.f90 const.o twobody_fun.o help.o define_types.o
	gfortran -c -O3 $< -fopenmp

twobody_fun.o : twobody_fun.f90 const.o help.o define_types.o distributions_fun.o
	gfortran -c -O3 $< -fopenmp

distributions_fun.o : distributions_fun.f90 const.o
	gfortran -c -O3 $< -fopenmp

data_out.o : data_out.f90 const.o
	gfortran -c $< -fopenmp

data_in.o : data_in.f90 const.o
	gfortran -c $< -fopenmp

help.o : help.f90 const.o
	gfortran -c -O3 $< -fopenmp

define_types.o : define_types.f90 const.o
	gfortran -c $<

const.o : const.f90 
	gfortran -c $< -fopenmp

clean :
	rm *.mod *.o *dudihc


