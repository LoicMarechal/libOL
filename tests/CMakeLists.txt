cmake_minimum_required(VERSION 3.9)

project(TEST_LIBOL)

include_directories( ${CMAKE_SOURCE_DIR}/../sources )
add_executable(basicTests basicTests.cpp ${CMAKE_SOURCE_DIR}/../sources )

# add extra checks
set(CMAKE_CXX_FLAGS_DEBUG "--std=c++11 -D_GLIBCXX_DEBUG -D_GLIBCXX_ASSERTIONS -fsanitize=address -fsanitize=undefined")

set(CMAKE_C_FLAGS_DEBUG "-D_GLIBCXX_DEBUG -D_GLIBCXX_ASSERTIONS -fsanitize=address -fsanitize=undefined")
