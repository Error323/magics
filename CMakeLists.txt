project (magics CXX C)
cmake_minimum_required (VERSION 2.8)

#-------------------------------------------------------------------------------
# General Settings
#-------------------------------------------------------------------------------
set (CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11")
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)

#-------------------------------------------------------------------------------
# Include sources
#-------------------------------------------------------------------------------
set (SOURCES_CPU
  src/timer.h
  src/bitboard.cpp
  src/island.cpp
  src/main.cpp
)

#-------------------------------------------------------------------------------
# Define executable and link libraries
#-------------------------------------------------------------------------------
add_executable (${CMAKE_PROJECT_NAME} ${SOURCES_CPU})
target_link_libraries (${CMAKE_PROJECT_NAME}
  rt
  pthread
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
message (STATUS "")
message (STATUS " CMakeflags (${CMAKE_PROJECT_NAME}):")
message (STATUS "   CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
message (STATUS "")
