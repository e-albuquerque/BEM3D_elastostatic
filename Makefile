FC = gfortran
FFLAGS = -O3 -fopenmp -fdefault-real-8
LIBS = -llapack -lblas
TARGET = quadesth_modern.exe

OBJS = data_mod.o solver_mod.o assembly_mod.o stress_mod.o io_mod.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -static $(LIBS) -o $(TARGET)

data_mod.o: data_mod.f90
	$(FC) $(FFLAGS) -c data_mod.f90

solver_mod.o: solver_mod.f90 data_mod.o
	$(FC) $(FFLAGS) -c solver_mod.f90

assembly_mod.o: assembly_mod.f90 data_mod.o
	$(FC) $(FFLAGS) -c assembly_mod.f90

stress_mod.o: stress_mod.f90 data_mod.o
	$(FC) $(FFLAGS) -c stress_mod.f90

io_mod.o: io_mod.f90 data_mod.o assembly_mod.o
	$(FC) $(FFLAGS) -c io_mod.f90

main.o: main.f90 data_mod.o solver_mod.o assembly_mod.o stress_mod.o io_mod.o
	$(FC) $(FFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(TARGET)
