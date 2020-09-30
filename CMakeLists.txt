cmake_minimum_required(VERSION 3.12)
project(raytracer)

# path needs to be up to where dso is placed (no recursive search)
find_library( EMBREE_LIB embree3 HINT $ENV{HOME}/toolchains/*/*)
find_path(EMBREE_INCLUDE_DIR embree3/rtcore.h HINT $ENV{HOME}/toolchains/*/*)

message(STATUS ${EMBREE_LIB})
message(STATUS $ENV{HOME}/toolchains/embree)
message(STATUS ${EMBREE_INCLUDE_DIR})

add_executable(a.out main.cpp)

target_include_directories(a.out BEFORE PUBLIC ${EMBREE_INCLUDE_DIR})
# target_link_directories(a.out BEFORE PUBLIC $ENV{HOME/toolchains/embree/lib})



target_link_libraries(a.out ${EMBREE_LIB})
