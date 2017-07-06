# 
# -------------  Makefile for Gemini & Dispec & Totido  ---------------
#

COMPILER=gcc
COMPILE_FLAGS= -c -O

#COMPILER=pgcc
#COMPILE_FLAGS= -c -O -mpentium

#################### Targets ###################
#
all: gemini dispec totido
#
#
gemini: 
	cd Gemini; make CC=$(COMPILER) CFLAGS="$(COMPILE_FLAGS)"; mv gemini ..
#
dispec:
	cd Dispec; make CC=$(COMPILER) CFLAGS="$(COMPILE_FLAGS)"; mv dispec ..
#
totido: 
	cd Totido; make CC=$(COMPILER) CFLAGS="$(COMPILE_FLAGS)"; mv totido ..
#
clean:
	-rm gemini totido dispec
	cd Gemini; rm -f *.o
	cd Dispec; rm -f *.o
	cd Totido; rm -f *.o
#
#
	-rm -f $(output)/d.*.sff
	echo \
	'set stalist=`genstationlist.tcl -e $(epimax) -gg2gc $(info)/ehb/info.$(EV).LHZ` \n \
	foreach s ($$stalist) \n \
	 	set datafile=`find $(datadir)/$(EV) -name \*$$s.\*.LHZ.sff -maxdepth 1` \n \
		echo "$$s $$datafile[1]" \n \
	 	set locopt=`locidopt.tcl $$datafile[1] $$s LHZ` \n \
		raw2procdata -hp $(highpass) -lp $(lowpass) -fu 5. -fo 100. -ct 100 $$locopt \\\n \
		   -s -i -r $(respdir) -tlen 1800 -dec $(ndec) -gg2gc \\\n \
		   -o $(output)/d.$$s.sff $$datafile[1] $(info)/ehb/info.$(EV).LHZ \n \
	end' | /bin/tcsh
	gensec.tcl -o $@ $(info)/ehb/info.$(EV).LHZ $(output)/d sff

