# Exact Diag
#FC     = /opt/mpich1/bin/mpif90
#FC     = ifort
FC     = mpif90
#FC     = mpifrtpx
##FC     = openmpif90
OPT    = -Wa -O3 -parallel -xSSE4.2 -zero -i-static -ip -ipo -ftz -openmp
#OPT    = -Wa -O3 -parallel -xCORE-AVX2 -zero -ip -ftz -qopenmp
#OPT    = -Wa -OO -parallel -ip -ipo -zero -i-static -openmp -traceback -check
#OPT    = -Wall -O3 -fopenmp 
#OPT    = -Wall -O0 -g -fopenmp
#OPT    = -Wa -O3 -parallel -xCORE-AVX2 -zero -ip -ftz -qopenmp
#OPT    = -Kfast,parallel,openmp,simd=2 -O3
#OPT    = -O3 -Kns,mfunc,eval,prefetch_conditional,ilfunc,fp_contract,fp_relaxed,openmp,parallel,simd=2 
RCHECK = -fno-range-check
OPT2   = -O0 -traceback -check 
OPT3   = -O3 -i-static -ip -ipo -ftz -zero -xSSE4.2
PROG   = Hybrid_TPS_
OBJS   = params.o output.o rand.o basic.o mpilib.o parallel.o Lanczos.o cg.o dsf.o main.o omplib.o timer.o check.o TPS.o

a.out: $(OBJS)
	$(FC) $(OPT) $(RCHECK) -o $@ $(OBJS)
	rm -f *mod $(OBJS)
main.o: $(PROG)main.f90 params.o timer.o rand.o basic.o Lanczos.o cg.o mpilib.o parallel.o dsf.o timer.o omplib.o output.o check.o TPS.o
	$(FC) -c $(OPT) -o $@ $<
TPS.o: $(PROG)TPS.f90 params.o rand.o basic.o mpilib.o parallel.o timer.o omplib.o output.o
	$(FC) -c $(OPT) -o $@ $<
dsf.o: $(PROG)DSF.f90 params.o mpilib.o parallel.o timer.o omplib.o output.o
	$(FC) -c $(OPT) -o $@ $<
cg.o: $(PROG)cg.f90 params.o rand.o basic.o mpilib.o parallel.o timer.o omplib.o output.o
	$(FC) -c $(OPT) -o $@ $<
Lanczos.o: $(PROG)Lanczos.f90 params.o rand.o basic.o mpilib.o parallel.o timer.o omplib.o output.o
	$(FC) -c $(OPT) -o $@ $<
output.o: $(PROG)output.f90 mpilib.o params.o
	$(FC) -c $(OPT) -o $@ $<
check.o: $(PROG)check.f90 params.o rand.o mpilib.o timer.o omplib.o parallel.o
	$(FC) -c $(OPT) -o $@ $<
parallel.o: $(PROG)parallel.f90 params.o rand.o mpilib.o timer.o omplib.o
	$(FC) -c $(OPT) -o $@ $<
mpilib.o: $(PROG)MPIlib.f90 params.o rand.o timer.o omplib.o
	$(FC) -c $(OPT) -o $@ $<

basic.o: $(PROG)basic.f90
	$(FC) -c $(OPT) -o $@ $<
params.o: $(PROG)params.f90
	$(FC) -c $(OPT) -o $@ $<
rand.o: $(PROG)ran.f90
	$(FC) $(RCHECK) -c $(OPT) -o $@ $<
omplib.o: $(PROG)OMPlib.f90
	$(FC) $(RCHECK) -c $(OPT) -o $@ $<
timer.o: $(PROG)timer.f90
	$(FC) -c $(OPT) -o $@ $<

remove :
	rm -f Salp*tp basisme* vec*tp seedme*
eval : $(PROG)eval.f90
	$(FC) $(OPT2) -o $@ $<
	rm -f *mod
eval_mpi : $(PROG)eval_mpi.f90
	$(FC) $(OPT3) -o $@ $<
	rm -f *mod
clean: 
	rm -f a.out $(OBJS) *mod list_fort eval
