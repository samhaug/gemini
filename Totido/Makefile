###################################################################### 
# -------------  Makefile for Totido  ---------------
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
TOTIDO = totido.o sft.o hbutcw.o writeAscii.o
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
#  You may omit this library if you comment out the call to writesff
#  in totido.f (close to the end)
#
#SFFLIB=-L. -lf2cstuff
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
totido: $(TOTIDO)
	$(COMPILER) -o $@ $(TOTIDO) $(SFFLIB) $(COMPLIB)
#
clean:
	-rm *.o gemini totido dispec
