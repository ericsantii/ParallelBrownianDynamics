#COMP= gfortran
#COMP= ifort -static-libcxa
#COMP= ifort -openmp
COMP= ifort -g -traceback -qopenmp
#COMP= ifort 
#  for intel f90 for P III
#GFLAGS= -O3 -xK -tpp6 -pad -Zp8 -wp_ipo 
#GFLAGS= -O3 -xK -tpp6 -pad -Zp8 -wp_ipo
#GFLAG77= -O3 -xK -tpp6 -pad -Zp8 -wp_ipo
#GFLAGSmain= -O2 -xK -tpp6 -pad -Zp8 -wp_ipo
#
#  for intel f90 for P IV  

#GFLAGS= -O3 -g -xW -tpp7 -pad -Zp8

#GFLAGS = -O3 -ipo -static
brow.x: vec.o part.o ll.o cell.o glob.o vfun.o quat.o matrices.o ran.o init.o gen.o sph.o dyn.o clus.o col.o out.o main.o 
	
	$(COMP) $(GFLAGS) vec.o part.o ll.o cell.o glob.o vfun.o quat.o matrices.o ran.o init.o gen.o sph.o dyn.o clus.o col.o out.o main.o -o brow.x

vec.o: vec.f90
	$(COMP) $(GFLAGS) vec.f90 -c

part.o: part.f90
	$(COMP) $(GFLAGS) part.f90 -c

ll.o: ll.f90
	$(COMP) $(GFLAGS) ll.f90 -c

cell.o: cell.f90
	$(COMP) $(GFLAGS) cell.f90 -c

glob.o: glob.f90
	$(COMP) $(GFLAGS) glob.f90 -c

vfun.o: vfun.f90
	$(COMP) $(GFLAGS) vfun.f90 -c

quat.o: quat.f90
	$(COMP) $(GFLAGS) quat.f90 -c

matrices.o: matrices.f90
	$(COMP) $(GFLAGS) matrices.f90 -c

ran.o: ran.f90
	$(COMP) $(GFLAGS) ran.f90 -c

init.o: init.f90
	$(COMP) $(GFLAGS) init.f90 -c

gen.o: gen.f90
	$(COMP) $(GFLAGS) gen.f90 -c

sph.o: sph.f90
	$(COMP) $(GFLAGS) sph.f90 -c

dyn.o: dyn.f90
	$(COMP) $(GFLAGS) dyn.f90 -c

clus.o: clus.f90
	$(COMP) $(GFLAGS) clus.f90 -c

col.o: col.f90
	$(COMP) $(GFLAGS) col.f90 -c

out.o: out.f90
	$(COMP) $(GFLAGS) out.f90 -c

main.o: main.f90
	$(COMP) $(GFLAGS) main.f90 -c

clean:
	rm -f *.o
	rm -f *.mod





