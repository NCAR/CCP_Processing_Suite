HOST = $(shell hostname)
null =
FC   = $(null)
#
# 'lens' at ORNL
ifeq ($(HOST),lens0)
  FC     = pgf95
  FFLAGS = -c
  NCDIRI = -I/ccs/home/wgstrand/include
  NCINCF = -lnetcdf
  NCDIRL = -L/ccs/home/wgstrand/lib
  NCLIBS = $(NCDIRL) $(NCINCF)
endif
#
# 'mineral' at NCAR/CGD
ifeq ($(HOST),mineral)
  FC     = pgf95
  FFLAGS = -c
  NCDIRI = -I/usr/local/include
  NCINCF = -lnetcdf
  NCDIRL = -L/usr/local/lib
  NCLIBS = $(NCDIRL) $(NCINCF)
endif
#
OBJCOMN = \
   array_tools.o \
   char_case.o \
   const_mod.o \
   kinds_mod.o \
   msg_mod.o \
   nf_tools.o \
   nf_wrap.o \
   nf_wrap_stubs.o \
   vars_mod.o

OBJCOCN = main_ccsm_ocn.o ccsm_ocn_rgd.o

OBJCICE = main_ccsm_ice.o ccsm_ice_rgd.o

all: ccsm_ice_rgd ccsm_ocn_rgd

ccsm_ice_rgd: $(OBJCOMN) $(OBJCICE)
	$(FC) -o ccsm_ice_rgd.$(HOST) $(OBJCOMN) $(OBJCICE) $(NCLIBS)

ccsm_ocn_rgd: $(OBJCOMN) $(OBJCOCN)
	$(FC) -o ccsm_ocn_rgd.$(HOST) $(OBJCOMN) $(OBJCOCN) $(NCLIBS)

clean:
	rm -f *.o *.mod

.SUFFIXES:            # Delete the default suffixes
.SUFFIXES: .F90 .o    # Define our suffix list

.F90.o:
	$(FC) $(FFLAGS) $(NCDIRI) $<

main_ccsm_ocn.o: main_ccsm_ocn.F90 kinds_mod.o vars_mod.o ccsm_ocn_rgd.o

main_ccsm_ice.o: main_ccsm_ice.F90 kinds_mod.o vars_mod.o ccsm_ice_rgd.o

ccsm_ocn_rgd.o: ccsm_ocn_rgd.F90 kinds_mod.o vars_mod.o msg_mod.o char_case.o const_mod.o nf_wrap.o nf_tools.o

ccsm_ice_rgd.o: ccsm_ice_rgd.F90 kinds_mod.o vars_mod.o msg_mod.o char_case.o const_mod.o nf_wrap.o nf_tools.o

array_tools.o: array_tools.F90 kinds_mod.o msg_mod.o

char_case.o: char_case.F90

const_mod.o: const_mod.F90 kinds_mod.o

kinds_mod.o: kinds_mod.F90

msg_mod.o: msg_mod.F90 kinds_mod.o

nf_tools.o: nf_tools.F90 kinds_mod.o array_tools.o const_mod.o nf_wrap.o

nf_wrap.o: nf_wrap.F90 kinds_mod.o msg_mod.o nf_wrap_stubs.o

nf_wrap_stubs.o: nf_wrap_stubs.F90 kinds_mod.o

vars_mod.o: vars_mod.F90 kinds_mod.o
