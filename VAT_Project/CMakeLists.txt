cmake_minimum_required (VERSION 2.6)
project (VAT)

find_package(Boost
  1.53.0 REQUIRED
  COMPONENTS program_options timer iostreams thread
  system timer chrono)

set(CMAKE_BUILD_TYPE Release)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -g -O3 -flto -march=native -fopenmp -Wno-narrowing -Wno-uninitialized -Wunused-variable -Wsign-compare -Wreorder -Wpsabi -Wsign-compare")
add_subdirectory("src/algo/blast/core")
include_directories(${CMAKE_SOURCE_DIR}/src/tools/include/avx256)
# include_directories(
#   "${CMAKE_SOURCE_DIR}/src"
#   "${Boost_INCLUDE_DIR} ")

add_executable(VAT src/VATMain.cpp
  src/basic/VATParameters.cpp
  src/tools/tinythread.cpp
  src/tools/sort_simd/common.cpp
  src/tools/sort_simd/avx256/merge_util.cpp
  src/tools/sort_simd/avx256/simd_sort.cpp
  src/tools/sort_simd/avx256/sort_util.cpp
  src/tools/sort_simd/avx256/Utils.cpp
  )

target_link_libraries(VAT blast_core ${Boost_LIBRARIES})

install(TARGETS VAT DESTINATION bin)
