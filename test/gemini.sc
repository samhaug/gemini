#!/bin/sh 
#######################################################################
#
#         Example-Script for GEMINI-program
#
#               J.R. Dalkolmo, April 1997
#######################################################################

Redirect_to=/dev/tty
#Redirect_to=gemini.sc.out

#-------------------------------------------------------------------
#                Set program parameters
#                ~~~~~~~~~~~~~~~~~~~~~~

# Choose the kind of motion you like to calculate. Set the
# variable to 1 for P-SV-motion only; set it to 2 for SH-motion only;
# set it to 3 to calculate all P-SV- and SH-motion.
What_Motion=$11

# Verbose level for monitor-output. '0' yields monitoring every 
# frequency; n every degree if mod(l,n)=0
Print_Level=0

# Length of seismogram in seconds
Seismo_Length=$3

# Damping time for complex frequency (->Laplace transform). A good
# choice is a fifth of the seismogram length.
Damping_Time=$4

# Minimum frequency in millihertz. You can set it to zero, 'Gemini' will
# adjust it to at least 1/Seismo_Length
Minimum_Frequency=$5

# Maximum frequency in millihertz
Maximum_Frequency=$6

# Take into account dispersion (attenuation), which leads to 
# frequency-dependent elastic moduli. Set 1 for 'yes', 0 for 'no'.
Dispersion_Switch=$7

# Minimum degree of spherical harmonics. We recommend '0'.
Minimum_Degree=$8

# Maximum degree of spherical harmonics. 'Gemini' will not compute
# beyond this limit.
Maximum_Degree=$9

# Step in the degree-domain. Normally this is '1', because of the
# 2*Pi-periodicity. Setting this greater than '1' speeds up the
# calculation and gives nice alias-effects because the Earth will
# become 2*Pi/lstep-periodic!!
Degree_Step=1

# Depth of the source in [km]. You can set this to Zero, if you like
Source_Depth=$10

# Accuracy which rules the performance of the integration 
# algorithm (Bulirsch-Stoer). Don't be too greedy, 1.e-4 should be
# sufficient
Accuracy=1.e-4

# File name of the earth model.
Earth_Model=$1

# File name of the window in the frequency-degree-domain, which contains
# tabulated maximum degrees for selected frequencies.
Omega_Ell_Window=$2

# Name of output file with expansion coefficients
Output_Filename=green/bas.f50.d100.3.out

# Confirm the input
Confirmation=1
#-------------------------------------------------------------------

#                            Run GEMINI
#                            ~~~~~~~~~~

if [ $# -a "$1" = -n ]; then
 /bin/rm $Output_Filename
 shift
fi

../gemini << HERE  >| $Redirect_to
$What_Motion
$Print_Level
$Seismo_Length
$Damping_Time
$Minimum_Frequency
$Maximum_Frequency
$Dispersion_Switch
$Minimum_Degree
$Maximum_Degree
$Degree_Step
$Source_Depth
$Accuracy
$Earth_Model
$Omega_Ell_Window
$Output_Filename
$Confirmation
HERE


