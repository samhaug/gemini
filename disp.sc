#!/bin/sh 
#######################################################################
#
#         Example-Script for DISPEC-program
#
#               J.R. Dalkolmo, August 1997
#######################################################################

#-------------------------------------------------------------------
#                Set program parameters
#                ~~~~~~~~~~~~~~~~~~~~~~


# Give here the file with the basis solutions calculated by GEMINI
BasisSolutions=green/bas.f50.d100.3.out

# The following variable must tell the source file containing ONE 
# set of earthquake parameters in Harvard-CMT-format
SourceFile=$1

# This defines the source mechanism: moment tensor ('m') 
# or single force ('f'):
SourceMech=m

# Set maximum order 'm' in the sum over spherical harmonics here.
# For example, set 'm=0' if the source is explosion type
MaximumOrder=$4

# This file must contain the window in the omega-l-domain also needed
# by GEMINI
OmEllWindow=$3

# Give here the length of the taper which is applied on the l-range
# to avoid cut-off-effects. This is an empirical number, don't choose
# it too large, you may cut in the surface wave branch
EllTaper=40

# This is the name of the file with station names and parameters defined
# in IRIS-DMC format
StationFile=$2

# Specify here          1 -> Station by latitude and longitude in degrees
# the receiver type:    2 -> Station by its abbreviation (e.g. PFO)
#                       3 -> all stations in file
#                       4 -> section along great circle between two points on sphere
ReceiverType=3

# Now, coordinates or station name or all stations or section?
# Example:     Type 1 :  Receiver='48.3 8.3'
#              Type 2 :  Receiver='BFO'
#              Type 3 :  Receiver='all'  (this is a dummy, 
#                                         program reads format '(1x)'
#              Type 4 :  Receiver='-30. -71. 48.3 8.3 10'
#
Receiver='all'

# Horizontal displacement is calculated by default in source centered
# coordiantes. If you like to have station centered coordinates, i.e.
# North-South- and East-West-coordinates give here a number greater
# than zero.
NSEWCoordinates=1
#
#  The output file will have the name "spec3k" in the current directory.
#
#-----------------------------------------

../dispec  << ENDE 
$BasisSolutions
$SourceFile
$SourceMech
$MaximumOrder
$OmEllWindow
$EllTaper
$StationFile
$ReceiverType
$Receiver
$NSEWCoordinates
ENDE
