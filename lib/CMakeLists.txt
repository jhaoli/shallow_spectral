cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(sources ./adquad.f ./fft99f.f ./ncarg.f ./netcdf.f)

add_library(fftpack ${sources})
