cmake_minimum_required(VERSION 3.22)

project(HighPerformancePricingEngine
        VERSION 1.0.0
        DESCRIPTION "Tier-1 High-Performance Derivatives Pricing Framework"
        LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Global C++ Configuration
# ------------------------------------------------------------------------------
# Strictly enforce modern C++20 standard without compiler-specific extensions
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

# Export compile commands for Clangd IntelliSense and VS Code integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Build Configurations & Compiler Flags
# ------------------------------------------------------------------------------
# Default to RelWithDebInfo (Optimized with debug symbols) if unspecified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
endif()

# Tier-1 hardware optimizations for low-latency quantitative computing
add_compile_options(
    -Wall -Wextra -Wpedantic -Werror # Treat all warnings as compilation errors
    $<$<CONFIG:Release>:-O3 -march=native -ffast-math -flto>
    $<$<CONFIG:RelWithDebInfo>:-O2 -march=native -g>
)

# ------------------------------------------------------------------------------
# Core Pricing Library Target
# ------------------------------------------------------------------------------
# Gather all source files from the src directory (Heston.cpp, pde_solvers.cpp, etc.)
file(GLOB_RECURSE CORE_SOURCES "src/*.cpp")

# Compile the core math and pricing engines as a static library
add_library(pricing_core_lib STATIC ${CORE_SOURCES})

# Publicly expose the include directory to anything linking against this target
target_include_directories(pricing_core_lib 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Link standard math and threading libraries required for Monte Carlo execution
find_package(Threads REQUIRED)
target_link_libraries(pricing_core_lib PRIVATE Threads::Threads m)

# ------------------------------------------------------------------------------
# Subdirectories (Bindings, Tests, Benchmarks)
# ------------------------------------------------------------------------------

# Python Pybind11 integration layer
option(BUILD_BINDINGS "Build Python Pybind11 dynamic library extension" ON)
if(BUILD_BINDINGS)
    add_subdirectory(bindings)
endif()

# Google Test Suite for mathematical validation
option(BUILD_TESTING "Build the C++ unit testing tree" ON)
if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()

# Google Benchmark Suite for latency regression tracking
option(BUILD_BENCHMARKS "Build the performance benchmarking tree" ON)
if(BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()