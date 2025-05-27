cmake_minimum_required(VERSION 3.16)
project(liborbit LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add library sources
file(GLOB LIB_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/Sources/*.cpp"
)
file(GLOB LIB_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Sources/*.h"
)

add_library(liborbit STATIC ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(liborbit PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Sources
)

# Optionally, add bindings subdirectory if present
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/bindings/CMakeLists.txt")
    add_subdirectory(bindings)
endif()

# Optionally, add executables
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Executables")
    file(GLOB EXEC_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Executables/*.cpp")
    foreach(EXEC_FILE ${EXEC_SOURCES})
        get_filename_component(EXEC_NAME ${EXEC_FILE} NAME_WE)
        add_executable(${EXEC_NAME} ${EXEC_FILE})
        target_link_libraries(${EXEC_NAME} PRIVATE liborbit)
        target_include_directories(${EXEC_NAME} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/Sources
        )
    endforeach()
endif()
