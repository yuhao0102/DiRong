SHELL           = /bin/sh
FC = mpiifort


# SOURCE FILES
#
SRCS_F90 = parse_namelist_mod.f90 spmd_init_mod.f90 grid_init_mod.f90 decomp_init_mod.f90 variable_init_mod.f90 coupling_ocn_model_mod.f90 model_setting_mod.f90 atm_demo.f90 
# SRCS_F90    = $(wildcard ./*.f90)
OBJS_ALL    = $(SRCS_F90:.f90=.o)

MPEUPATH = /public3/home/sc52823/mct_test_model/model_platform/models/libs/MCT/mpeu
MCTPATH = /public3/home/sc52823/mct_test_model/model_platform/models/libs/MCT/mct

#MPIPATH = /opt/intel/impi/3.2.0.011
#NETCDFPATH = /opt/netCDF-intel13-without-hdf5
MPIPATH=/public1/soft/intel/2018/compilers_and_libraries_2018.2.199/linux/mpi/intel64/
NETCDFPATH=/public1/soft/netcdf/4.4.1-parallel-icc18/


MCTLIBS    = -isystem $(MPIPATH)/include64 -I$(NETCDFPATH)/include -L$(NETCDFPATH)/lib -lnetcdff -lnetcdf -I$(MPEUPATH) -L$(MPEUPATH) -I$(MCTPATH) -L$(MCTPATH) -lmct -lmpeu

all: atm_demo

atm_demo: $(OBJS_ALL)
	$(FC) -o $@ $(OBJS_ALL) $(FCFLAGS) $(UTLDFLAGS) $(MCTLIBS) $()

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	    $(FC) -c -o $@ $(MCTLIBS) $<

clean:
	    ${RM} *.o *.mod atm_demo

