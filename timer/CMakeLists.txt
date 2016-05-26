
cmake_minimum_required(VERSION 3.0)
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${Boost_INCLUDE_DIRS}")

if(Boost_FOUND)
  add_library(timer timer.cpp)
  target_link_libraries(timer ${Boost_LIBRARIES})
  message(STATUS "CMakeLists in timer")
endif()
