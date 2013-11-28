project (magics CXX C)
cmake_minimum_required (VERSION 2.8)

#-------------------------------------------------------------------------------
# General Settings
#-------------------------------------------------------------------------------
set (CMAKE_CXX_FLAGS "-Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

option (ENABLE_CUDA "Enable/disable cuda" ON)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)

if (ENABLE_CUDA)
  find_package (CUDA)
endif (ENABLE_CUDA)

#-------------------------------------------------------------------------------
# Include sources
#-------------------------------------------------------------------------------
set (SOURCES_CPU
  src/timer.h
  src/ga.cpp
  src/random.c
)

set (SOURCES_GPU
  src/timer.h
  src/Verbose.cpp
  src/main.cpp
  src/cuda/init.h
  src/cuda/init.cu
  src/cuda/process.h
  src/cuda/process.cu
)

#-------------------------------------------------------------------------------
# Define executable and link libraries
#-------------------------------------------------------------------------------
if (CUDA_FOUND)
  include_directories("/usr/local/cuda/include")
  link_directories("/usr/local/cuda/lib64")
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -gencode arch=compute_30,code=sm_30 -Xptxas=-v)
  cuda_add_executable (${CMAKE_PROJECT_NAME}-gpu ${SOURCES_GPU})
  target_link_libraries (${CMAKE_PROJECT_NAME}-gpu
    rt
    curand
  )
endif (CUDA_FOUND)

add_executable (${CMAKE_PROJECT_NAME}-cpu ${SOURCES_CPU})
target_link_libraries (${CMAKE_PROJECT_NAME}-cpu
  rt
)

#-------------------------------------------------------------------------------
# Status report
#-------------------------------------------------------------------------------
message (STATUS "")
message (STATUS " C/C++:")
message (STATUS "   C++ Compiler:                ${CMAKE_CXX_COMPILER}")
message (STATUS "   C++  flags (Release):        ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
message (STATUS "   C++  flags (Debug):          ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
message (STATUS "   C++  flags (RelWithDebInfo): ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message (STATUS "   NVCC flags:                  ${CUDA_NVCC_FLAGS}")
message (STATUS "")
message (STATUS " CMakeflags (${CMAKE_PROJECT_NAME}):")
message (STATUS "   CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
message (STATUS "   ENABLE_CUDA                 ${ENABLE_CUDA}")
message (STATUS "")
