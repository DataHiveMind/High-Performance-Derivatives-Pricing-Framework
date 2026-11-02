# tests/CMakeLists.txt

# Fetch Google Test if not already provided by the system
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

# Define the testing executable
add_executable(test_models 
    test_models.cpp
    test_greeks.cpp
)

# Link against the core pricing library and gtest
target_link_libraries(test_models PRIVATE 
    pricing_core_lib 
    gtest_main 
    gmock 
    Threads::Threads
)

# Enable CTest integration and automatic test discovery
include(GoogleTest)
gtest_discover_tests(test_models)