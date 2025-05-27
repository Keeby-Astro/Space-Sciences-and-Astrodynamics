cmake_minimum_required(VERSION 3.16)

# Add all C++ test sources in this directory
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Use Catch2 if available
find_package(Catch2 3 QUIET)

foreach(TEST_FILE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_include_directories(${TEST_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/curtis_scripts/cpp
        ${CMAKE_SOURCE_DIR}/curtis_scripts/cpp/Sources
    )
    if(Catch2_FOUND)
        target_link_libraries(${TEST_NAME} PRIVATE Catch2::Catch2WithMain)
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endif()
endforeach()
