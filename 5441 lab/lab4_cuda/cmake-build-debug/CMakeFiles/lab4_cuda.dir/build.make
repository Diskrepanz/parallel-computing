# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/shengding/Desktop/5441 lab/lab4_cuda"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/lab4_cuda.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lab4_cuda.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lab4_cuda.dir/flags.make

CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.o: CMakeFiles/lab4_cuda.dir/flags.make
CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.o: ../Sheng_Ding_cuda.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.o"
	/usr/local/Cellar/gcc/9.2.0_1/bin/gcc-9 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.o   -c "/Users/shengding/Desktop/5441 lab/lab4_cuda/Sheng_Ding_cuda.c"

CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.i"
	/usr/local/Cellar/gcc/9.2.0_1/bin/gcc-9 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/Users/shengding/Desktop/5441 lab/lab4_cuda/Sheng_Ding_cuda.c" > CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.i

CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.s"
	/usr/local/Cellar/gcc/9.2.0_1/bin/gcc-9 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/Users/shengding/Desktop/5441 lab/lab4_cuda/Sheng_Ding_cuda.c" -o CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.s

# Object files for target lab4_cuda
lab4_cuda_OBJECTS = \
"CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.o"

# External object files for target lab4_cuda
lab4_cuda_EXTERNAL_OBJECTS =

lab4_cuda: CMakeFiles/lab4_cuda.dir/Sheng_Ding_cuda.c.o
lab4_cuda: CMakeFiles/lab4_cuda.dir/build.make
lab4_cuda: CMakeFiles/lab4_cuda.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable lab4_cuda"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lab4_cuda.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lab4_cuda.dir/build: lab4_cuda

.PHONY : CMakeFiles/lab4_cuda.dir/build

CMakeFiles/lab4_cuda.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lab4_cuda.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lab4_cuda.dir/clean

CMakeFiles/lab4_cuda.dir/depend:
	cd "/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/shengding/Desktop/5441 lab/lab4_cuda" "/Users/shengding/Desktop/5441 lab/lab4_cuda" "/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug" "/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug" "/Users/shengding/Desktop/5441 lab/lab4_cuda/cmake-build-debug/CMakeFiles/lab4_cuda.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/lab4_cuda.dir/depend
