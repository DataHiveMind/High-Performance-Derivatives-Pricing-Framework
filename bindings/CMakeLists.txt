# bindings/CMakeLists.txt

# Fetch Pybind11 using FetchContent for seamless dependency management
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# Define the Python module target
pybind11_add_module(pricing_engine python_bindings.cpp)

# Link the core C++ pricing library to the Python module
target_link_libraries(pricing_engine PRIVATE pricing_core_lib)

# Apply extreme optimization flags for the Python bridge
target_compile_options(pricing_engine PRIVATE -O3 -flto -march=native)