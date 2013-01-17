#
# Makefile for CESM CMOR codeb40.abrupt4xco2_solar.beta19.001_cam2.h0_process.sh
#
MACH = $(shell echo $(HOSTNAME) | cut -c1-6)
null =
FC   = $(null)
#
# 'silver' at NCAR/CGD
#
ifeq ($(MACH),silver)
  NCFDIR  = /usr/local/netcdf-4.1.3-pgi-hpf-cc-11.5-0
  UDUDIR  = /home/strandwg/proj/software/udunits-2.1.24
  UIDDIR  = /datalocal/ccpa/strandwg/software/uuid-1.6.2
  CMRDIR  = /datalocal/ccpa/strandwg/CMIP5/CMOR2/current
  HD5DIR  = /usr/local/hdf5-1.8.7
  KRBDIR  = /usr/kerberos/lib64
#
  NCFINC  = -I${NCFDIR}/include -I./include
  NCFLIB  = -L${NCFDIR}/lib -lnetcdf -lnetcdff
  KRBLIB  = -L${KRBDIR}/lib -lcurl -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lidn -lz
  UDUINC  = -I${UDUDIR}/include
  UDULIB  = -L${UDUDIR}/lib -ludunits2 -lexpat
  UIDINC  = -I${UIDDIR}/include
  UIDLIB  = -L${UIDDIR}/lib -luuid
  FC	  = /usr/local/pgi/linux86/bin/pgf95
#  FFLAGS  = -V -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -Mbounds -fpic
  FFLAGS  = -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -Mbounds -fpic
endif
#
# 'modeling1' at NCAR/ACD
#
ifeq ($(MACH),modeli)
  NCFDIR  = /usr/local/netcdf-4.1.2
  UDUDIR  = /data1/strandwg/software/udunits-2.1.21
  UIDDIR  = /data1/strandwg/software/uuid-1.6.2
  CMRDIR  = /data1/strandwg/software/cmor
  HD5DIR  = /usr/lib64
  KRBDIR  = /usr/lib64
#
  NCFINC  = -I${NCFDIR}/include -I./include
  NCFLIB  = -L${NCFDIR}/lib -lnetcdf -lnetcdff
  KRBLIB  = -L${KRBDIR} -lcurl -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lidn -lssl -lcrypto -lz
  UDUINC  = -I${UDUDIR}/include
  UDULIB  = -L${UDUDIR}/lib -ludunits2 -lexpat
  UIDINC  = -I${UIDDIR}/include
  UIDLIB  = -L${UIDDIR}/lib -luuid
  FC      = /usr/pgi/linux86-64/11.3/bin/pgf95
  FFLAGS  = -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -Mbounds -fpic
#
  HD5INC  = -I${HD5DIR}/include
  HD5LIB  = -L${HD5DIR} -lhdf5_hl -lhdf5 -lz -lm
endif
#
# 'geyser' CISL
#
ifeq ($(MACH),yslogi)
  NCFDIR  = /glade/apps/opt/netcdf/4.2/pgi/default
  UDUDIR  = /glade/apps/opt/udunits/2.1.24/pgi/12.5
  UIDDIR  = /glade/u/home/strandwg/software/uuid-1.6.2
  CMRDIR  = /glade/home/strandwg/software/cmor
  HD5DIR  = /glade/apps/opt/hdf5/1.8.9/pgi/default
  KRBDIR  = /usr/lib64
#
  NCFINC  = -I${NCFDIR}/include -I./include
  NCFLIB  = -L${NCFDIR}/lib -lnetcdf -lnetcdff
  KRBLIB  = -L${KRBDIR}/lib -lcurl -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lidn -lssl -lcrypto -lz
  UDUINC  = -I${UDUDIR}/include
  UDULIB  = -L${UDUDIR}/lib -ludunits2 -lexpat
  UIDINC  = -I${UIDDIR}/include
  UIDLIB  = -L${UIDDIR}/lib -luuid
  FC	  = pgf95
  FFLAGS  = -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -fpic
endif
#
# 'euclid' at NERSC
#
ifeq ($(MACH),euclid)
  NCFDIR  = /usr/common/usg/netcdf/4.1.1/pgi
  UDUDIR  = /usr/common/usg/udunits/default/pgi
  UIDDIR  = /global/u2/s/strandwg/software/uuid-1.6.2
  CMRDIR  = /global/u2/s/strandwg/software/cmor
  HD5DIR  = /usr/common/usg/hdf5/1.8.6/pgi/serial
  KRBDIR  = /usr/lib64
#
  NCFINC  = -I${NCFDIR}/include -I./include
  NCFLIB  = -L${NCFDIR}/lib -lnetcdf -lnetcdff
#  KRBLIB  = -L${KRBDIR} -lcurl -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lidn -lssl -lcrypto -lz
  KRBLIB  = -L${KRBDIR} -ldl -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lssl -lcrypto -lz
  UDUINC  = -I${UDUDIR}/include
  UDULIB  = -L${UDUDIR}/lib -ludunits2 -lexpat
  UIDINC  = -I${UIDDIR}/include
  UIDLIB  = -L${UIDDIR}/lib -luuid
  FC	  = /usr/common/usg/pgi/12.3/linux86-64/12.3/bin/pgf95
#  FFLAGS  = -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -Mbounds -m64 -mcmodel medium
  FFLAGS  = -Wl,-rpath=${UDUDIR}/lib -Wl,-rpath=${UIDDIR}/lib -Mbounds -fpic
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
	  define_lnd_axes.o \
	  define_atm_axes.o \
	  define_ice_axes.o \
	  define_ocn_axes.o \
	  load_expt.o \
	  reset_netcdf_var.o \
	  load_xwalk.o \
	  add_global_metadata.o \
	  file_output_times.o \
	  parse_rip.o \
	  dphybrid.o \
	  vertintp.o \
	  parse_ncfile.o \
	  get_ice_grid.o \
	  get_ocn_grid.o \
	  get_lnd_grid.o \
	  get_atm_grid.o

all:	  Amon_CMOR 3hrSlev_CMOR 3hrMlev_CMOR 3hrCurt_CMOR 3hrPlev_CMOR sites_CMOR day_CMOR Lmon_CMOR \
	  LImon_CMOR OImon_CMOR Omon_CMOR fx_CMOR 6hrLev_CMOR 6hrPlev_CMOR 3hr_CMOR.f90 aero_CMOR \
	  cfMon_CMOR cfDay_CMOR cf3hr_CMOR Oyr_CMOR cfSites_CMOR

Amon_CMOR: $(OBJ) Amon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Amon_CMOR.o \
	-o $(EXEDIR)/Amon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

cfMon_CMOR: $(OBJ) cfMon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) cfMon_CMOR.o \
	-o $(EXEDIR)/cfMon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

cfDay_CMOR: $(OBJ) cfDay_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) cfDay_CMOR.o \
	-o $(EXEDIR)/cfDay_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

cf3hr_CMOR: $(OBJ) cf3hr_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) cf3hr_CMOR.o \
	-o $(EXEDIR)/cf3hr_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

aero_CMOR: $(OBJ) aero_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) aero_CMOR.o \
	-o $(EXEDIR)/aero_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

3hrSlev_CMOR: $(OBJ) 3hrSlev_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 3hrSlev_CMOR.o \
	-o $(EXEDIR)/3hrSlev_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

3hrMlev_CMOR: $(OBJ) 3hrMlev_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 3hrMlev_CMOR.o \
	-o $(EXEDIR)/3hrMlev_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

3hrCurt_CMOR: $(OBJ) 3hrCurt_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 3hrCurt_CMOR.o \
	-o $(EXEDIR)/3hrCurt_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

sites_CMOR: $(OBJ) sites_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) sites_CMOR.o \
	-o $(EXEDIR)/sites_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

cfSites_CMOR: $(OBJ) cfSites_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) cfSites_CMOR.o \
	-o $(EXEDIR)/cfSites_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

3hrPlev_CMOR: $(OBJ) 3hrPlev_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 3hrPlev_CMOR.o \
	-o $(EXEDIR)/3hrPlev_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

6hrLev_CMOR: $(OBJ) 6hrLev_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 6hrLev_CMOR.o \
	-o $(EXEDIR)/6hrLev_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

6hrPlev_CMOR: $(OBJ) 6hrPlev_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 6hrPlev_CMOR.o \
	-o $(EXEDIR)/6hrPlev_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

3hr_CMOR: $(OBJ) 3hr_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) 3hr_CMOR.o \
	-o $(EXEDIR)/3hr_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

fx_CMOR: $(OBJ) fx_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) fx_CMOR.o \
	-o $(EXEDIR)/fx_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

Lmon_CMOR: $(OBJ) Lmon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Lmon_CMOR.o \
	-o $(EXEDIR)/Lmon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

LImon_CMOR: $(OBJ) LImon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) LImon_CMOR.o \
	-o $(EXEDIR)/LImon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

OImon_CMOR: $(OBJ) OImon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) OImon_CMOR.o \
	-o $(EXEDIR)/OImon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

Omon_CMOR: $(OBJ) Omon_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Omon_CMOR.o \
	-o $(EXEDIR)/Omon_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

Oyr_CMOR: $(OBJ) Oyr_CMOR.o
	$(FC) $(FFLAGS) $(NCFINC) $(HD5INC) $(UDUINC) $(UIDINC) $(OBJ) Oyr_CMOR.o \
	-o $(EXEDIR)/Oyr_CMOR $(UDULIB) $(NCFLIB) $(UIDLIB) $(CMRLIB) $(HD5LIB) $(KRBLIB) $(UDULIB)

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

dphybrid.o: dphybrid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c dphybrid.f90

get_atm_grid.o: get_atm_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_atm_grid.f90

get_ice_grid.o: get_ice_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_ice_grid.f90

get_ocn_grid.o: get_ocn_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_ocn_grid.f90

get_lnd_grid.o: get_lnd_grid.f90 mymods.o
	$(FC) $(FFLAGS) $(NCFINC) -c get_lnd_grid.f90

build_filenames.o: build_filenames.f90 mymods.o
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

reset_netcdf_var.o: reset_netcdf_var.f90 module_netcdf.o netcdf_jfl.o
	$(FC) $(FFLAGS) $(LOCINC) -c reset_netcdf_var.f90

parse_rip.o: parse_rip.f90
	$(FC) $(FFLAGS) $(LOCINC) -c parse_rip.f90

parse_ncfile.o: parse_ncfile.f90
	$(FC) $(FFLAGS) $(LOCINC) -c parse_ncfile.f90

module_netcdf.o: module_netcdf.f90
	$(FC) $(FFLAGS) $(NCFINC) -c module_netcdf.f90

netcdf_jfl.o: netcdf_jfl.f90 module_netcdf.o
	$(FC) $(FFLAGS) $(NCFINC) -c netcdf_jfl.f90

handle_err.o: handle_err.f90
	$(FC) $(FFLAGS) $(NCFINC) -c handle_err.f90

Amon_CMOR.o: Amon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Amon_CMOR.f90

cfMon_CMOR.o: cfMon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c cfMon_CMOR.f90

cfDay_CMOR.o: cfDay_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c cfDay_CMOR.f90

cf3hr_CMOR.o: cf3hr_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c cf3hr_CMOR.f90

aero_CMOR.o: aero_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c aero_CMOR.f90

3hrSlev_CMOR.o: 3hrSlev_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 3hrSlev_CMOR.f90

3hrMlev_CMOR.o: 3hrMlev_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 3hrMlev_CMOR.f90

3hrCurt_CMOR.o: 3hrCurt_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 3hrCurt_CMOR.f90

sites_CMOR.o: sites_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c sites_CMOR.f90

cfSites_CMOR.o: cfSites_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c cfSites_CMOR.f90

3hrPlev_CMOR.o: 3hrPlev_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 3hrPlev_CMOR.f90

6hrLev_CMOR.o: 6hrLev_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 6hrLev_CMOR.f90

6hrPlev_CMOR.o: 6hrPlev_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 6hrPlev_CMOR.f90

3hr_CMOR.o: 3hr_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c 3hr_CMOR.f90

fx_CMOR.o: fx_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c fx_CMOR.f90

Lmon_CMOR.o: Lmon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Lmon_CMOR.f90

LImon_CMOR.o: LImon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c LImon_CMOR.f90

OImon_CMOR.o: OImon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c OImon_CMOR.f90

Omon_CMOR.o: Omon_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Omon_CMOR.f90

Oyr_CMOR.o: Oyr_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c Oyr_CMOR.f90

day_CMOR.o: day_CMOR.f90 $(OBJ)
	$(FC) $(FFLAGS) $(NCFINC) $(LOCINC) $(CMRINC) $(UIDINC) -c day_CMOR.f90

clean: 
	/bin/rm -f *.o *.mod

