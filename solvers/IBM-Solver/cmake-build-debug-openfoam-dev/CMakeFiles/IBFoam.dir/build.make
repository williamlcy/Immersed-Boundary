# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/william/clion-2023.3.2/bin/cmake/linux/x64/bin/cmake

# The command to remove a file.
RM = /home/william/clion-2023.3.2/bin/cmake/linux/x64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev

# Include any dependencies generated for this target.
include CMakeFiles/IBFoam.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/IBFoam.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/IBFoam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/IBFoam.dir/flags.make

CMakeFiles/IBFoam.dir/IBFoam.C.o: CMakeFiles/IBFoam.dir/flags.make
CMakeFiles/IBFoam.dir/IBFoam.C.o: /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/IBFoam.C
CMakeFiles/IBFoam.dir/IBFoam.C.o: CMakeFiles/IBFoam.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/IBFoam.dir/IBFoam.C.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/IBFoam.dir/IBFoam.C.o -MF CMakeFiles/IBFoam.dir/IBFoam.C.o.d -o CMakeFiles/IBFoam.dir/IBFoam.C.o -c /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/IBFoam.C

CMakeFiles/IBFoam.dir/IBFoam.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/IBFoam.dir/IBFoam.C.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/IBFoam.C > CMakeFiles/IBFoam.dir/IBFoam.C.i

CMakeFiles/IBFoam.dir/IBFoam.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/IBFoam.dir/IBFoam.C.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/IBFoam.C -o CMakeFiles/IBFoam.dir/IBFoam.C.s

# Object files for target IBFoam
IBFoam_OBJECTS = \
"CMakeFiles/IBFoam.dir/IBFoam.C.o"

# External object files for target IBFoam
IBFoam_EXTERNAL_OBJECTS =

IBFoam: CMakeFiles/IBFoam.dir/IBFoam.C.o
IBFoam: CMakeFiles/IBFoam.dir/build.make
IBFoam: CMakeFiles/IBFoam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable IBFoam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/IBFoam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/IBFoam.dir/build: IBFoam
.PHONY : CMakeFiles/IBFoam.dir/build

CMakeFiles/IBFoam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/IBFoam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/IBFoam.dir/clean

CMakeFiles/IBFoam.dir/depend:
	cd /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev /home/william/OpenFOAM/william-v2312/20240403_IBM/IBM-Solver/cmake-build-debug-openfoam-dev/CMakeFiles/IBFoam.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/IBFoam.dir/depend

