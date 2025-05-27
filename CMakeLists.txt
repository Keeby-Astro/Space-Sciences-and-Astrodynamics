cmake_minimum_required(VERSION 3.26)
project(space_mech LANGUAGES CXX)

add_subdirectory(curtis_scripts/cpp)
enable_testing()
add_subdirectory(tests/cpp)
