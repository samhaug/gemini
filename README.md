# gemini

This directory contains the Gemini code of Joerg Dalkolmo and Wolfgang Friederich.
I have made some minor alterations to make this code more user friendly.

setup\_bin contains gemini\_station & ndk\_2\_gemini, which allow you to 
convert an NDK file to a gemini source file as well as make a gemini station
file from a pickled obspy stream.

The executable do\_all is a concise way to run all of the gemini code. It copies
all of the output files into one directory as well as some of the input files 
that are useful for remembering the details of a run.
