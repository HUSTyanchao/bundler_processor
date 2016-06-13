
cmake_minimum_required(VERSION 3.0)

project(bundler_processer)

if(WIN32)  
  message(STATUS "Now is windows")
  set("PROJECT_BINARY_DIR" "${PROJECT_SOURCE_DIR}/vsbuild")  
elseif(APPLE)  
  message(STATUS "Now is Apple systens.")  
elseif(UNIX)  
  message(STATUS "Now is UNIX-like OS's. Including APPLE OS X and CygWin") 
  set("PROJECT_BINARY_DIR" "${PROJECT_SOURCE_DIR}/cygbuild")  
endif()  

# The version number.
set(HELLO_VERSION_MAJOR 0)


# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")

# add opencv library
set(OpenCV_DIR "D:/Program Files/opencv/build")
find_package(OpenCV REQUIRED)
set(EXTRA_LIBS ${EXTRA_LIBS} ${OpenCV_LIBS})

# add boost library
set(Boost_INCLUDE_DIR "D:/boost_1_61_0")
set(Boost_LIBRARY_DIR "D:/boost_1_61_0/lib64-msvc-12.0")
find_package(Boost REQUIRED COMPONENTS Chrono Timer)
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
   link_directories(${Boost_LIBRARY_DIR})
   message(STATUS "boost found")
else()
   message(FATAL_ERROR "boost not found")
endif()


# source and header of the timer
set(timer_SRC timer/timer.cpp)
set(timer_HDR timer/timer.h)

option(USE_EXIF_READER "use exif_reader" ON)
# source and header of the exif reader
if(USE_EXIF_READER)
include_directories("${PROJECT_SOURCE_DIR}/exif_reader")
add_subdirectory(exif_reader)
set(EXTRA_LIBS ${EXTRA_LIBS} exif_reader)
endif()

# add the executable
aux_source_directory(. DIR_SRCS)
add_executable(bundler_processer ${DIR_SRCS} ${exif_SRC} ${exif_HDR} ${timer_SRC} ${timer_HDR})
target_link_libraries(bundler_processer ${EXTRA_LIBS} ${Boost_LIBRARIES})

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )