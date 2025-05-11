cmake_minimum_required(VERSION 3.10)

project(PAIN CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Collect all source files
file(GLOB SOURCES "*.cpp")
list(APPEND SOURCES "../tests/tests.cpp")

# Create the executable from all source files
add_executable(test.out ${SOURCES})
