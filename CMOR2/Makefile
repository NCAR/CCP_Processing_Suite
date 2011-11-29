#
# Makefile for CESM CMOR code
#
MACH = $(shell echo $(HOST) | cut -c1-6)
null =
FC   = $(null)
#
# 'copper' at NCAR/CGD
#
ifeq ($(MACH),copper)
  NCFDIR  = /usr/local/netcdf-4.1.1-pgi-hpf-cc-7.5-2
  UDUDIR  = /datalocal/ccpa/strandwg/software/udunits-2.1.21
  UIDDIR  = /datalocal/ccpa/strandwg/software/uuid-1.6.2
  CMRDIR  = /datalocal/ccpa/strandwg/CMIP5/CMOR2/current
  HD5DIR  = /usr/local/hdf5
  KRBDIR  = /usr/kerberos/lib64
#
  NCFINC  = -I${NCFDIR}/include -I./include
  NCFLIB  = -L${NCFDIR}/lib -lnetcdf
  KRBLIB  = -L${KRBDIR}/lib -lcurl -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lidn -lssl -lcrypto -lz
  UDUINC  = -I${UDUDIR}/include
  UDULIB  = -L${UDUDIR}/lib -ludunits2 -lexpat
  UIDINC  = -I${UIDDIR}/include
  UIDLIB  = -L${UIDDIR}/lib -luuid
  FC	  = /usr/local/pgi/linux86/bin/pgf95
  FFLAGS  = -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -module ${CMRDIR}
endif
#
# 'mirage' cluser in CISL
ifeq ($(MACH),mirage)
  NCFDIR  = /fs/local/apps/netcdf-4.1.3_intel
  UDUDIR  = /fs/local/apps/udunits-2.1.19
  UIDDIR  = /glade/home/strandwg/cmor2/uuid-1.6.2
  CMRDIR  = /glade/home/strandwg/cmor2/cmor
  HD5DIR  = /fs/local/apps/hdf5-1.8.7
  KRBDIR  = /usr/lib64
#
  NCFINC  = -I${NCFDIR}/include -I./include
  NCFLIB  = -L${NCFDIR}/lib -lnetcdf -lnetcdff
  KRBLIB  = -L${KRBDIR}/lib -lcurl -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lidn -lssl -lcrypto -lz
  UDUINC  = -I${UDUDIR}/include
  UDULIB  = -L${UDUDIR}/lib -ludunits2 -lexpat
  UIDINC  = -I${UIDDIR}/include
  UIDLIB  = -L${UIDDIR}/lib -luuid
  FC	  = /fs/local/bin/ifort
endif
#
#
CMRINC  = -I${CMRDIR}/include
CMRLIB  = -L${CMRDIR}/lib -lcmor
HD5INC  = -I${HD5DIR}/include
HD5LIB  = -L${HD5DIR}/lib -lhdf5_hl -lhdf5 -lz -lm
NCFINC  = -I${NCFDIR}/include -I./include
#
LOCINC  = -I.
#
SRCDIR  = $(HOME)/CCP_Processing_Suite/CMOR2
EXEDIR  = $(HOME)/bin
OBJ	= mymods.o \
	  handle_err.o \
	  module_netcdf.o \
	  netcdf_jfl.o \
	  get_cmor_info.o \
	  get_exp_metadata.o \
 	  build_filenames.o \
	  define_atm_axes.o \
	  define_lnd_axes.o \
	  define_ice_axes.o \
	  define_ocn_axes.o \
	  load_expt.o \
	  load_table_info.o \
	  load_xwalk.o \
	  add_global_metadata.o \
	  file_output_times.o \
	  parse_rip.o \
	  vertintp.o \
	  parse_ncfile.o \
	  get_ice_grid.o \
	  get_ocn_grid.o \
	  get_lnd_grid.o \
	  get_atm_grid.o

all:	  Amon_CMOR day_CMOR Lmon_CMOR OImon_CMOR Omon_CMOR fx_CMOR

Amon_CMOR: $(OBJ) Amon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Amon_CMOR.o \
	-o $(EXEDIR)/Amon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

fx_CMOR: $(OBJ) fx_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) fx_CMOR.o \
	-o $(EXEDIR)/fx_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

Lmon_CMOR: $(OBJ) Lmon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Lmon_CMOR.o \
	-o $(EXEDIR)/Lmon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

OImon_CMOR: $(OBJ) OImon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) OImon_CMOR.o \
	-o $(EXEDIR)/OImon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

Omon_CMOR: $(OBJ) Omon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Omon_CMOR.o \
	-o $(EXEDIR)/Omon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

day_CMOR: $(OBJ) day_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) day_CMOR.o \
	-o $(EXEDIR)/day_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

mymods.o: mymods.f90
	$(FC) $(FFLAGS) $(NCFINC) -c mymods.f90

get_cmor_info.o: get_cmor_info.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) -c get_cmor_info.f90

get_exp_metadata.o: get_exp_metadata.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) -c get_exp_metadata.f90

add_global_metadata.o: add_global_metadata.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) $(CMRINC) -c add_global_metadata.f90

load_expt.o: load_expt.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) -c load_expt.f90

file_output_times.o: file_output_times.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) -c file_output_times.f90

vertintp.o: vertintp.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c vertintp.f90

get_atm_grid.o: get_atm_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_atm_grid.f90

get_ice_grid.o: get_ice_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_ice_grid.f90

get_ocn_grid.o: get_ocn_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_ocn_grid.f90

get_lnd_grid.o: get_lnd_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_lnd_grid.f90

build_filenames.o: build_filenames.f90
	$(FC) $(FFLAGS) $(LOCINC) -c build_filenames.f90

define_atm_axes.o: define_atm_axes.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) $(CMRINC) -c define_atm_axes.f90

define_ice_axes.o: define_ice_axes.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) $(CMRINC) -c define_ice_axes.f90

define_ocn_axes.o: define_ocn_axes.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) $(CMRINC) -c define_ocn_axes.f90

define_lnd_axes.o: define_lnd_axes.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) $(CMRINC) -c define_lnd_axes.f90

load_xwalk.o: load_xwalk.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) -c load_xwalk.f90

load_table_info.o: load_table_info.f90 mymods.o
	$(FC) $(FFLAGS) $(LOCINC) -c load_table_info.f90

parse_rip.o: parse_rip.f90
	$(FC) $(FFLAGS) $(LOCINC) -c parse_rip.f90

parse_ncfile.o: parse_ncfile.f90
	$(FC) $(FFLAGS) $(LOCINC) -c parse_ncfile.f90

module_netcdf.o: module_netcdf.f90
	$(FC) $(FFLAGS) $(NCFINC) -c module_netcdf.f90

netcdf_jfl.o: netcdf_jfl.f90
	$(FC) $(FFLAGS) $(NCFINC) -c netcdf_jfl.f90

handle_err.o: handle_err.f90
	$(FC) $(FFLAGS) $(NCFINC) -c handle_err.f90

Amon_CMOR.o: Amon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Amon_CMOR.f90

fx_CMOR.o: fx_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c fx_CMOR.f90

Lmon_CMOR.o: Lmon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Lmon_CMOR.f90

OImon_CMOR.o: OImon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c OImon_CMOR.f90

Omon_CMOR.o: Omon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Omon_CMOR.f90

day_CMOR.o: day_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c day_CMOR.f90

clean: 
	/bin/rm -f *.o *.mod

