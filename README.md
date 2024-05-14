# **sofatodat.c**

SofatoDAT.c is a commandline tool, to convert file from the sofa standard for stored HRTF files, to the .DAT standard required for use
with the CSound HRTF opcodes.

## **Description**

An in-depth paragraph about your project and overview of use.

## **Getting Started**

## **Dependencies**

The file requires the following libraries to be installed

<stdlib.h>

<stdio.h>

<string.h>

<unistd.h>

<netcdf.h>  //netcdf Library used to access data in sofa file

<math.h>

<fftw3.h>   //Fastest fourier transform in the west library used to interpolate data

some standard, though some more specific

## **Installing**

The .c file can be run from the command line, or terminal on a mac
It will ask for the location of the .sofa file to be converted
It will ask for the location of the .dat file to save to
.DAT files will be created with a sample rate and L/R appended to the name you select, allowing the created file
to work in CSound

The program will output a log, during the run, to aid in identifying issues, and can be useful for following the processes occuring.


**Authors**

Thom McDonnell

**Version History**

1.0

