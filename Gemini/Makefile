# #####################################################################
# -------------------  Makefile for Gemini   ---------------
#
# Compiler and flags are handed over from the makefile in the parent
# directory.
######################################################################

#-------------------Options and library/-paths -----------------------
F2CFLAGS= -u -f
FFLAGS= -e -O -ffixed-line-length-none
CFLAGS= -c
CLIB = -L/usr/lib -lf2c -lm 
FC=g77
#
# ------------------ Dependent o-files ------------------------------
#
GEMINI = gemini.o propag.o ixone.o aasmst.o gtllim.o \
	rzextr.o gtmanis.o stavani.o zetl.o cpanis.o odemin.o oderad.o\
	odepco.o odetor.o odeliq.o kodein.o bsstep.o mmid.o lyrfnd.o \
   stripmn.o rescal.o
#
#
#
#**************    Rules ************************************
#
COMPILER=g77
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
#
gemini: $(GEMINI)
	$(COMPILER) -o $@ $(GEMINI) $(COMPLIB)
#
clean:
	-rm *.o gemini
