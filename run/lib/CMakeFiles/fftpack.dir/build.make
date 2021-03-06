# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/qjz/work/study/shallow_spectral

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qjz/work/study/shallow_spectral/run

# Include any dependencies generated for this target.
include lib/CMakeFiles/fftpack.dir/depend.make

# Include the progress variables for this target.
include lib/CMakeFiles/fftpack.dir/progress.make

# Include the compile flags for this target's objects.
include lib/CMakeFiles/fftpack.dir/flags.make

lib/CMakeFiles/fftpack.dir/adquad.f.o: lib/CMakeFiles/fftpack.dir/flags.make
lib/CMakeFiles/fftpack.dir/adquad.f.o: ../lib/adquad.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qjz/work/study/shallow_spectral/run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object lib/CMakeFiles/fftpack.dir/adquad.f.o"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qjz/work/study/shallow_spectral/lib/adquad.f -o CMakeFiles/fftpack.dir/adquad.f.o

lib/CMakeFiles/fftpack.dir/adquad.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/fftpack.dir/adquad.f.i"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qjz/work/study/shallow_spectral/lib/adquad.f > CMakeFiles/fftpack.dir/adquad.f.i

lib/CMakeFiles/fftpack.dir/adquad.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/fftpack.dir/adquad.f.s"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qjz/work/study/shallow_spectral/lib/adquad.f -o CMakeFiles/fftpack.dir/adquad.f.s

lib/CMakeFiles/fftpack.dir/adquad.f.o.requires:

.PHONY : lib/CMakeFiles/fftpack.dir/adquad.f.o.requires

lib/CMakeFiles/fftpack.dir/adquad.f.o.provides: lib/CMakeFiles/fftpack.dir/adquad.f.o.requires
	$(MAKE) -f lib/CMakeFiles/fftpack.dir/build.make lib/CMakeFiles/fftpack.dir/adquad.f.o.provides.build
.PHONY : lib/CMakeFiles/fftpack.dir/adquad.f.o.provides

lib/CMakeFiles/fftpack.dir/adquad.f.o.provides.build: lib/CMakeFiles/fftpack.dir/adquad.f.o


lib/CMakeFiles/fftpack.dir/fft99f.f.o: lib/CMakeFiles/fftpack.dir/flags.make
lib/CMakeFiles/fftpack.dir/fft99f.f.o: ../lib/fft99f.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qjz/work/study/shallow_spectral/run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object lib/CMakeFiles/fftpack.dir/fft99f.f.o"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qjz/work/study/shallow_spectral/lib/fft99f.f -o CMakeFiles/fftpack.dir/fft99f.f.o

lib/CMakeFiles/fftpack.dir/fft99f.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/fftpack.dir/fft99f.f.i"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qjz/work/study/shallow_spectral/lib/fft99f.f > CMakeFiles/fftpack.dir/fft99f.f.i

lib/CMakeFiles/fftpack.dir/fft99f.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/fftpack.dir/fft99f.f.s"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qjz/work/study/shallow_spectral/lib/fft99f.f -o CMakeFiles/fftpack.dir/fft99f.f.s

lib/CMakeFiles/fftpack.dir/fft99f.f.o.requires:

.PHONY : lib/CMakeFiles/fftpack.dir/fft99f.f.o.requires

lib/CMakeFiles/fftpack.dir/fft99f.f.o.provides: lib/CMakeFiles/fftpack.dir/fft99f.f.o.requires
	$(MAKE) -f lib/CMakeFiles/fftpack.dir/build.make lib/CMakeFiles/fftpack.dir/fft99f.f.o.provides.build
.PHONY : lib/CMakeFiles/fftpack.dir/fft99f.f.o.provides

lib/CMakeFiles/fftpack.dir/fft99f.f.o.provides.build: lib/CMakeFiles/fftpack.dir/fft99f.f.o


lib/CMakeFiles/fftpack.dir/ncarg.f.o: lib/CMakeFiles/fftpack.dir/flags.make
lib/CMakeFiles/fftpack.dir/ncarg.f.o: ../lib/ncarg.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qjz/work/study/shallow_spectral/run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object lib/CMakeFiles/fftpack.dir/ncarg.f.o"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qjz/work/study/shallow_spectral/lib/ncarg.f -o CMakeFiles/fftpack.dir/ncarg.f.o

lib/CMakeFiles/fftpack.dir/ncarg.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/fftpack.dir/ncarg.f.i"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qjz/work/study/shallow_spectral/lib/ncarg.f > CMakeFiles/fftpack.dir/ncarg.f.i

lib/CMakeFiles/fftpack.dir/ncarg.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/fftpack.dir/ncarg.f.s"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qjz/work/study/shallow_spectral/lib/ncarg.f -o CMakeFiles/fftpack.dir/ncarg.f.s

lib/CMakeFiles/fftpack.dir/ncarg.f.o.requires:

.PHONY : lib/CMakeFiles/fftpack.dir/ncarg.f.o.requires

lib/CMakeFiles/fftpack.dir/ncarg.f.o.provides: lib/CMakeFiles/fftpack.dir/ncarg.f.o.requires
	$(MAKE) -f lib/CMakeFiles/fftpack.dir/build.make lib/CMakeFiles/fftpack.dir/ncarg.f.o.provides.build
.PHONY : lib/CMakeFiles/fftpack.dir/ncarg.f.o.provides

lib/CMakeFiles/fftpack.dir/ncarg.f.o.provides.build: lib/CMakeFiles/fftpack.dir/ncarg.f.o


lib/CMakeFiles/fftpack.dir/netcdf.f.o: lib/CMakeFiles/fftpack.dir/flags.make
lib/CMakeFiles/fftpack.dir/netcdf.f.o: ../lib/netcdf.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qjz/work/study/shallow_spectral/run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object lib/CMakeFiles/fftpack.dir/netcdf.f.o"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qjz/work/study/shallow_spectral/lib/netcdf.f -o CMakeFiles/fftpack.dir/netcdf.f.o

lib/CMakeFiles/fftpack.dir/netcdf.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/fftpack.dir/netcdf.f.i"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qjz/work/study/shallow_spectral/lib/netcdf.f > CMakeFiles/fftpack.dir/netcdf.f.i

lib/CMakeFiles/fftpack.dir/netcdf.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/fftpack.dir/netcdf.f.s"
	cd /home/qjz/work/study/shallow_spectral/run/lib && /usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qjz/work/study/shallow_spectral/lib/netcdf.f -o CMakeFiles/fftpack.dir/netcdf.f.s

lib/CMakeFiles/fftpack.dir/netcdf.f.o.requires:

.PHONY : lib/CMakeFiles/fftpack.dir/netcdf.f.o.requires

lib/CMakeFiles/fftpack.dir/netcdf.f.o.provides: lib/CMakeFiles/fftpack.dir/netcdf.f.o.requires
	$(MAKE) -f lib/CMakeFiles/fftpack.dir/build.make lib/CMakeFiles/fftpack.dir/netcdf.f.o.provides.build
.PHONY : lib/CMakeFiles/fftpack.dir/netcdf.f.o.provides

lib/CMakeFiles/fftpack.dir/netcdf.f.o.provides.build: lib/CMakeFiles/fftpack.dir/netcdf.f.o


# Object files for target fftpack
fftpack_OBJECTS = \
"CMakeFiles/fftpack.dir/adquad.f.o" \
"CMakeFiles/fftpack.dir/fft99f.f.o" \
"CMakeFiles/fftpack.dir/ncarg.f.o" \
"CMakeFiles/fftpack.dir/netcdf.f.o"

# External object files for target fftpack
fftpack_EXTERNAL_OBJECTS =

lib/libfftpack.a: lib/CMakeFiles/fftpack.dir/adquad.f.o
lib/libfftpack.a: lib/CMakeFiles/fftpack.dir/fft99f.f.o
lib/libfftpack.a: lib/CMakeFiles/fftpack.dir/ncarg.f.o
lib/libfftpack.a: lib/CMakeFiles/fftpack.dir/netcdf.f.o
lib/libfftpack.a: lib/CMakeFiles/fftpack.dir/build.make
lib/libfftpack.a: lib/CMakeFiles/fftpack.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qjz/work/study/shallow_spectral/run/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking Fortran static library libfftpack.a"
	cd /home/qjz/work/study/shallow_spectral/run/lib && $(CMAKE_COMMAND) -P CMakeFiles/fftpack.dir/cmake_clean_target.cmake
	cd /home/qjz/work/study/shallow_spectral/run/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fftpack.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/CMakeFiles/fftpack.dir/build: lib/libfftpack.a

.PHONY : lib/CMakeFiles/fftpack.dir/build

lib/CMakeFiles/fftpack.dir/requires: lib/CMakeFiles/fftpack.dir/adquad.f.o.requires
lib/CMakeFiles/fftpack.dir/requires: lib/CMakeFiles/fftpack.dir/fft99f.f.o.requires
lib/CMakeFiles/fftpack.dir/requires: lib/CMakeFiles/fftpack.dir/ncarg.f.o.requires
lib/CMakeFiles/fftpack.dir/requires: lib/CMakeFiles/fftpack.dir/netcdf.f.o.requires

.PHONY : lib/CMakeFiles/fftpack.dir/requires

lib/CMakeFiles/fftpack.dir/clean:
	cd /home/qjz/work/study/shallow_spectral/run/lib && $(CMAKE_COMMAND) -P CMakeFiles/fftpack.dir/cmake_clean.cmake
.PHONY : lib/CMakeFiles/fftpack.dir/clean

lib/CMakeFiles/fftpack.dir/depend:
	cd /home/qjz/work/study/shallow_spectral/run && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qjz/work/study/shallow_spectral /home/qjz/work/study/shallow_spectral/lib /home/qjz/work/study/shallow_spectral/run /home/qjz/work/study/shallow_spectral/run/lib /home/qjz/work/study/shallow_spectral/run/lib/CMakeFiles/fftpack.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/CMakeFiles/fftpack.dir/depend

