#
# Makefile para montar el programa DFT4He3d 
#  (Density Functional Theory 4He 3dimensional program)
#
# (Version de compilacion para Pentium IV)
#
COMP = ifort
CFLAGS = -c -fast -fpp -fopenmp -parallel -qopt-matmul -unroll 
#CFLAGS = -c -g -check all -static 
LD_FLAGS = -threads -parallel -qopt-matmul

#   Name of the program
PROGNAME=readwf-slurm

#   Fortran objects
objs=titols.o	readwf-slurm.o
#
.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(objs)
	$(COMP)	-o ../$(PROGNAME) $(objs)  $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;

clean:
	rm -f *.o *.bak *.lst *.mod;
distclean:
	rm -f *.o *.bak *.lst *.mod ../$(PROGNAME);
