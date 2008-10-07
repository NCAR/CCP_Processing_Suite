FC     = pgf95
FFLAGS = -c
NCDIRI = -I/usr/local/netcdf-3.6.1-pgi-hpf-cc-6.1-6/include
NCINCF = -lnetcdf
NCDIRL = -L/usr/local/netcdf-3.6.1-pgi-hpf-cc-6.1-6/lib
NCLIBS = $(NCDIRL) $(NCINCF)

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

OBJPOCN = main_pcm1_ocn.o pcm1_ocn_rgd.o

OBJCOCN = main_ccsm_ocn.o ccsm_ocn_rgd.o

OBJPICE = main_pcm1_ice.o pcm1_ice_rgd.o

OBJCICE = main_ccsm_ice.o ccsm_ice_rgd.o

all: pcm1_ice_rgd pcm1_ocn_rgd pcm1_ocn_trans ccsm_ice_rgd ccsm_ocn_rgd

pcm1_ice_rgd: $(OBJCOMN) $(OBJPICE)
	$(FC) -o pcm1_ice_rgd $(OBJCOMN) $(OBJPICE) $(NCLIBS)

pcm1_ocn_rgd: $(OBJCOMN) $(OBJPOCN)
	$(FC) -o pcm1_ocn_rgd $(OBJCOMN) $(OBJPOCN) $(NCLIBS)

pcm1_ocn_trans: pcm1_ocn_trans.F90
	$(FC) -o pcm1_ocn_trans $(NCDIRI) $(NCDIRL) pcm1_ocn_trans.F90 $(NCINCF)

ccsm_ice_rgd: $(OBJCOMN) $(OBJCICE)
	$(FC) -o ccsm_ice_rgd $(OBJCOMN) $(OBJCICE) $(NCLIBS)

ccsm_ocn_rgd: $(OBJCOMN) $(OBJCOCN)
	$(FC) -o ccsm_ocn_rgd $(OBJCOMN) $(OBJCOCN) $(NCLIBS)

clean:
	rm -f *.o *.mod

.SUFFIXES:            # Delete the default suffixes
.SUFFIXES: .F90 .o    # Define our suffix list

.F90.o:
	$(FC) $(FFLAGS) $(NCDIRI) $<

main_ccsm_ocn.o: main_ccsm_ocn.F90 kinds_mod.o vars_mod.o ccsm_ocn_rgd.o

main_ccsm_ice.o: main_ccsm_ice.F90 kinds_mod.o vars_mod.o ccsm_ice_rgd.o

main_pcm1_ocn.o: main_pcm1_ocn.F90 kinds_mod.o vars_mod.o pcm1_ocn_rgd.o

main_pcm1_ice.o: main_pcm1_ice.F90 kinds_mod.o vars_mod.o pcm1_ice_rgd.o

ccsm_ocn_rgd.o: ccsm_ocn_rgd.F90 kinds_mod.o vars_mod.o msg_mod.o char_case.o const_mod.o nf_wrap.o nf_tools.o

ccsm_ice_rgd.o: ccsm_ice_rgd.F90 kinds_mod.o vars_mod.o msg_mod.o char_case.o const_mod.o nf_wrap.o nf_tools.o

pcm1_ocn_rgd.o: pcm1_ocn_rgd.F90 kinds_mod.o vars_mod.o msg_mod.o char_case.o const_mod.o nf_wrap.o nf_tools.o

pcm1_ice_rgd.o: pcm1_ice_rgd.F90 kinds_mod.o vars_mod.o msg_mod.o char_case.o const_mod.o nf_wrap.o nf_tools.o

array_tools.o: array_tools.F90 kinds_mod.o msg_mod.o

char_case.o: char_case.F90

const_mod.o: const_mod.F90 kinds_mod.o

kinds_mod.o: kinds_mod.F90

msg_mod.o: msg_mod.F90 kinds_mod.o

nf_tools.o: nf_tools.F90 kinds_mod.o array_tools.o const_mod.o nf_wrap.o

nf_wrap.o: nf_wrap.F90 kinds_mod.o msg_mod.o nf_wrap_stubs.o

nf_wrap_stubs.o: nf_wrap_stubs.F90 kinds_mod.o

vars_mod.o: vars_mod.F90 kinds_mod.o
