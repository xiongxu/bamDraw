# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

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
CMAKE_COMMAND = /share/work1/staff/xuxiong/software/backup/cmake-3.0.0/bin/cmake

# The command to remove a file.
RM = /share/work1/staff/xuxiong/software/backup/cmake-3.0.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build

# Include any dependencies generated for this target.
include CMakeFiles/bam_draw_0_0_5_cmake.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bam_draw_0_0_5_cmake.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bam_draw_0_0_5_cmake.dir/flags.make

CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o: CMakeFiles/bam_draw_0_0_5_cmake.dir/flags.make
CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o: bam_draw_0.0.5.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o   -c /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build/bam_draw_0.0.5.c

CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build/bam_draw_0.0.5.c > CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.i

CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build/bam_draw_0.0.5.c -o CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.s

CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.requires:
.PHONY : CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.requires

CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.provides: CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.requires
	$(MAKE) -f CMakeFiles/bam_draw_0_0_5_cmake.dir/build.make CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.provides.build
.PHONY : CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.provides

CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.provides.build: CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o

# Object files for target bam_draw_0_0_5_cmake
bam_draw_0_0_5_cmake_OBJECTS = \
"CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o"

# External object files for target bam_draw_0_0_5_cmake
bam_draw_0_0_5_cmake_EXTERNAL_OBJECTS =

bam_draw_0_0_5_cmake: CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o
bam_draw_0_0_5_cmake: CMakeFiles/bam_draw_0_0_5_cmake.dir/build.make
bam_draw_0_0_5_cmake: CMakeFiles/bam_draw_0_0_5_cmake.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bam_draw_0_0_5_cmake"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bam_draw_0_0_5_cmake.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bam_draw_0_0_5_cmake.dir/build: bam_draw_0_0_5_cmake
.PHONY : CMakeFiles/bam_draw_0_0_5_cmake.dir/build

CMakeFiles/bam_draw_0_0_5_cmake.dir/requires: CMakeFiles/bam_draw_0_0_5_cmake.dir/bam_draw_0.0.5.c.o.requires
.PHONY : CMakeFiles/bam_draw_0_0_5_cmake.dir/requires

CMakeFiles/bam_draw_0_0_5_cmake.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bam_draw_0_0_5_cmake.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bam_draw_0_0_5_cmake.dir/clean

CMakeFiles/bam_draw_0_0_5_cmake.dir/depend:
	cd /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build /home/xuxiong/work/c/bam_handle/bam_draw_0_0_5_build/CMakeFiles/bam_draw_0_0_5_cmake.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bam_draw_0_0_5_cmake.dir/depend

