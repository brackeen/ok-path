# mkdir -p build && cd build && cmake -G Xcode ..
cmake_minimum_required(VERSION 3.6.0)

project(ok_path_test C)

include_directories(..)
file(GLOB ok_path_test_files "main.c" "../ok_path.c" "../ok_path.h")
add_executable(ok_path_test ${ok_path_test_files})
source_group("" FILES ${ok_path_test_files})

# Enable warnings
if (CMAKE_C_COMPILER_ID MATCHES "Clang")
  # Enable -Wwrite-strings because -Weverything doesn't enable it in all versions of Clang
  set_target_properties(ok_path_test PROPERTIES COMPILE_FLAGS "-Weverything -Wwrite-strings -Wno-padded")
elseif (CMAKE_C_COMPILER_ID MATCHES "GNU")
  # Disable unused-functions because of this GCC bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=64079
  set_target_properties(ok_path_test PROPERTIES COMPILE_FLAGS "-Wall -Wextra -Wwrite-strings -Wno-unused-function")
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
  # Disable 'struct padding'
  set_target_properties(ok_path_test PROPERTIES COMPILE_FLAGS "/Wall /wd4820 /wd4324")
endif()

# CTest setup, using valgrind if found.
enable_testing()
if (NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
  find_program(MEMCHECK_COMMAND valgrind)
  set(MEMCHECK_COMMAND_OPTIONS "--error-exitcode=1" "--leak-check=full")
endif()
if (NOT MEMCHECK_COMMAND)
  add_test(NAME ok_path_test COMMAND ok_path_test)
else()
  add_test(NAME ok_path_test_memcheck COMMAND ${MEMCHECK_COMMAND} ${MEMCHECK_COMMAND_OPTIONS} ./ok_path_test)
endif()
