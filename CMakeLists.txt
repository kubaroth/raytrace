cmake_minimum_required(VERSION 3.12)
project(raytracer)

# path needs to be up to where dso is placed (no recursive search)
find_library( EMBREE_LIB embree3 HINT $ENV{HOME}/toolchains/*/*)
find_path(EMBREE_INCLUDE_DIR embree3/rtcore.h HINT $ENV{HOME}/toolchains/*/*)

find_library( TBB_LIB tbb HINT $ENV{HOME}/toolchains/*/*)
find_path(TBB_INCLUDE_DIR tbb HINT $ENV{HOME}/toolchains/*/*)


message(STATUS ${EMBREE_LIB})
message(STATUS ${EMBREE_INCLUDE_DIR})

message(STATUS ${TBB_LIB})

message(STATUS ${TBB_INCLUDE_DIR})

message (STATUS ${EMBREE_DIR})

# Use custom dso location in installed target
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(a.out main.cpp)

target_include_directories(a.out BEFORE PUBLIC
  ${EMBREE_INCLUDE_DIR}
  # /home/kuba/SRC/embree/common/
  # /home/kuba/SRC/embree/common/math
  # /home/kuba/SRC/embree/common/sys
  # /home/kuba/SRC/embree/common/tasking
  ${TBB_INCLUDE_DIR}
  )

# to avoid problems with fininding tbb library at runinm (without explicitly
# setting LD_LIBRARY_PATH include header #include <tbb/tbb.h>
# This will mark library as NEED (readelf -d ./a.out) and will
# force linker to search specified directory
target_link_libraries(a.out  ${TBB_LIB} ${EMBREE_LIB})

install(TARGETS a.out)
