cmake_minimum_required(VERSION 3.20)

# Project name and version
project(ScientificComputingCode VERSION 0.1.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add the library
add_library(ScientificComputingLib
    ConjugateGradient.cpp
    GivensRotations.cpp
    Matrix.cpp
    Vec.cpp
)


# Add the main executable
add_executable(main main.cpp)
target_link_libraries(main PRIVATE ScientificComputingLib)

# Add Google Test
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Find all test files in the tests directory
file(GLOB TEST_SOURCES "tests/*.cpp")

# Add the test executable
add_executable(runTests ${TEST_SOURCES})
target_link_libraries(runTests PRIVATE ScientificComputingLib ${GTEST_LIBRARIES} pthread)

# Add a test target
add_test(NAME runTests COMMAND runTests)