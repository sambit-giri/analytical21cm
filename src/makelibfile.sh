#!/bin/sh
rm -f *.o *~ *.mod fort* *.nfs0* libANRAGfortlib.a ANRAGpylib.so
gfortran -c -O3 -ffast-math -fPIC -shared  nrtype.f90
gfortran -c -O3 -ffast-math -fPIC -shared  nrutil.f90
gfortran -c -O3 -ffast-math -fPIC -shared  nr.f90
gfortran -c -O3 -ffast-math -fPIC -shared  adaptint.f90
gfortran -c -O3 -ffast-math -fPIC -shared  rkck.f90
gfortran -c -O3 -ffast-math -fPIC -shared  rkqs.f90
gfortran -c -O3 -ffast-math -fPIC -shared  odeint.f90
gfortran -c -O3 -ffast-math -fPIC -shared  param.f90
gfortran -c -O3 -ffast-math -fPIC -shared  cosmo.f90
gfortran -c -O3 -ffast-math -fPIC -shared  press_sch.f90
gfortran -c -O3 -ffast-math -fPIC -shared  funct.f90
gfortran -c -O3 -ffast-math -fPIC -shared  subr.f90
gfortran -c -O3 -ffast-math -fPIC -shared  subr_FZH.f90
gfortran -c -O3 -ffast-math -fPIC -shared  subr_plot.f90
gfortran -c -O3 -ffast-math -fPIC -shared  subr_main.f90
gfortran -c -O3 -ffast-math -fPIC -shared  controler.f90
ar -rcs libANRAGfortlib.a nrtype.o nrutil.o nr.o adaptint.o rkck.o rkqs.o odeint.o param.o cosmo.o press_sch.o funct.o subr.o subr_FZH.o subr_plot.o subr_main.o
f2py -c controler.f90 -L. -lANRAGfortlib -m ANRAGpylib


