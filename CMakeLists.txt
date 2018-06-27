cmake_minimum_required(VERSION 3.7)
project(mx2)

set(CMAKE_CXX_STANDARD 11)

file(GLOB SOURCE_FILES
        "*.h"
        "*.cpp"
        "*.tpp"
)

add_executable(mx2.x ${SOURCE_FILES})
