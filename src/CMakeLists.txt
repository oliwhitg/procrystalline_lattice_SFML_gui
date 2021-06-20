cmake_minimum_required(VERSION 3.7)
project(PROCRYSTAL)

set(CMAKE_CXX_STANDARD 11)

##
set(SFML_INCLUDE_DIR "/home/SFML-2.5.1/include")
set(SFML_LIBRARY_DIR "/home/SFML-2.5.1/lib")
set(SFML_DIR "/home/SFML-2.5.1/lib/cmake/SFML")

link_directories(SFML_LIBRARY_DIR)
include_directories(SFML_INCLUDE_DIR)

find_package(SFML 2.5.1 COMPONENTS system window graphics network audio)

if(SFML_FOUND)
    message(STATUS "SFML_INCLUDE_DIR: ${SFML_INCLUDE_DIR}")
    message(STATUS "SFML_LIBRARIES: ${SFML_LIBRARIES}")
    message(STATUS "SFML_VERSION: ${SFML_VERSION}")
endif()
##

file(GLOB SOURCE_FILES
        "*.h"
        "lattice.cpp"
        "chain.cpp"
        "main.cpp"
        "node.cpp"
        "ring.cpp"
        "sampler.cpp"
        "*.tpp"
)

add_executable(procrystal.x ${SOURCE_FILES})

##
target_link_libraries (procrystal.x sfml-graphics sfml-window sfml-system)
##