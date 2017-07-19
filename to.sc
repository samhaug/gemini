#!/bin/sh
#
#                Script for TOTIDO
#
#              J.R. Dalkolmo, September 1997
#-----------------------------------------------------------------------
#if [ $# = 0 ]; then
#   echo usage: to.sc [-r \<responsefile\>] [-s \<timeshift\>] \
#   [-l \<lowpasscornerfrequency\>] [-L \<lowpassorder\>] \
#   [-h \<highpasscornerfrequency\>] [-H \<highpassorder\>] \
#   [-o \<secondsout\>] [-p \<zeropadding\>] [-O \<outputformat\>] \
#   -f spectrumfile
#   exit 1
#fi

# Give here the name of the file with the spectra generated by dispec

# Response file containing real and imaginary part of instrument
# transfer function at the frequencies used in gemini and dispec.
# Each line of the file contains freqeuncy, real part, imaginary part.
ResponseFile=""

#  Time by which beginning of times series is delayed
TimeShift='0.0'

# To apply a Butterworth-filter (low and high pass) 
# to the seismogram specify in the string
# the number n of filters and then n times  order and corner frequency.
# For example: One low pass with order 7 and corner frequency 0.01 Hz
#              Two high pass: order 3 and corner freq. 0.008 Hz AND 
#                             order 2 and corner freq. 0.007 Hz
#             -------> type: 1  7  0.01
#                            2  3  0.008  2 0.007
LowpassNumber=$4
LowpassOrder=$5
LowpassCF=$6
HighpassNumber=$7
HighpassOrder=$8
HighpassCF=$9

#  either a for Ascii output or s for SFF output
OutputFormat="a"

# The Fast Fourier Transform needs 2**n samples. If the spectra do not have
# such a length zeros will be appended to accomplish this. If you give here
# a number n greater than zero, the number of samples will be multiplicated 
# with 2**n, resulting in interpolation/smoothing of the time series.
ZeroPadding=2

# With the following character you can choose the type of seismogram:
#     type d for displacement, v for velocity, a for acceleration or
#     g for accelerometer response.
SeismoType=$3

# Specify here the length of the output time series in seconds 
SecondsOut=$2

#Specify the spectrum file here
SpectrumFile=$1

#opts="L:l:H:h:s:o:O:p:r:f:"

#while getopts $opts actopt;  
# do
#    case $actopt in
#         L) LowpassOrder=$OPTARG;;
#         l) LowpassCF=$OPTARG;;
#         H) HighpassOrder=$OPTARG;;
#         h) HighpassCF=$OPTARG;;
#         s) TimeShift=$OPTARG;;
#         o) SecondsOut=$OPTARG;;
#         O) OutputFormat=$OPTARG;;
#         p) ZeroPadding=$OPTARG;;
#         r) ResponseFile=$OPTARG;;
#         f) SpectrumFile=$OPTARG;;
#    esac
# done


#      Running Totido

$10/totido << END
$SpectrumFile
$ResponseFile
$TimeShift
$LowpassNumber $LowpassOrder $LowpassCF
$HighpassNumber $HighpassOrder $HighpassCF
$ZeroPadding
$SeismoType
$SecondsOut
$OutputFormat
END

#stuplo plot_z.sff
