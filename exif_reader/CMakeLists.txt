
cmake_minimum_required(VERSION 3.0)
include_directories("${PROJECT_BINARY_DIR}")
message("exif_reader: ${PROJECT_BINARY_DIR}")

set(exif_SRC exif_reader.cpp  jhead-2.90/exif.cpp jhead-2.90/gpsinfo.cpp jhead-2.90/iptc.cpp jhead-2.90/jhead.cpp jhead-2.90/jpgfile.cpp jhead-2.90/makernote.cpp jhead-2.90/paths.cpp )
set(exif_HDR exif_reader.h    jhead-2.90/jhead.h)

add_library(exif_reader ${exif_SRC} ${exif_HDR})
