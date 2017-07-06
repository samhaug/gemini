# 
######################################################################
# ---------------------  Makefile for Dispec  ---------------
#
# Compiler and flags are handed over from the makefile in the parent
# directory.
######################################################################
#
#-------------------Options and library/-paths -----------------------
F2CFLAGS= -u -f
FFLAGS= -e -O
CFLAGS= -c
CLIB = -L/usr/lib -lf2c -lm 
#
# ------------------ Dependent o-files ------------------------------
#
DISPEC = dispec.o epitra.o sharm.o  momfor.o gtllim.o \
	costap.o elapa.o stnorm.o getStation.o \
      readcmt.o readSForce.o necklace.o allStations.o
#
#**************    Rules ************************************
#
COMPILER=f77
#
# We prefer working with 'f2c & gcc', which generate smaller
# and faster executables. If you don't like this or you haven't
# got the programs above, just comment out the following lines.
# The library CLIB is only used by 'f2c & gcc'.
#
# Comment out from here ........................ 
COMPILER=$(CC)
COMPLIB=$(CLIB)
#
%.o: %.c
.f.o:
	@f2c $(F2CFLAGS) $<
	$(COMPILER) $(CFLAGS) $(CPPFLAGS) $(<:.f=.c) -o $@
	@rm  $(<:.f=.c)
# ...................................  to here.
#
#
#################### Targets ###################
#
dispec: $(DISPEC)
	$(COMPILER) -o $@ $(DISPEC) $(COMPLIB)
#
clean:
	-rm *.o gemini totido dispec
