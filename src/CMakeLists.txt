cmake_minimum_required(VERSION 3.30)

set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH}
    /usr/lib64/cmake/Catch2/)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Catch2)
find_library(GMP_LIBRARY gmp)

project(
    algebraic_complex_numbers
    VERSION 0.1
    LANGUAGES C CXX
)

add_library(algebraic_complex_numbers STATIC arith.cpp matrix.cpp)

add_executable(tests tests.cpp arith.cpp matrix.cpp)
target_link_libraries(tests PRIVATE algebraic_complex_numbers gmp)
