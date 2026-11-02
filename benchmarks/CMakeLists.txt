# benchmarks/CMakeLists.txt

find_package(benchmark REQUIRED)

# Monte Carlo Throughput Benchmark
add_executable(bench_monte_carlo bench_monte_carlo.cpp)

# Link against the core pricing library, Google Benchmark, and threading (pthread)
target_link_libraries(bench_monte_carlo PRIVATE 
    pricing_core_lib 
    benchmark::benchmark 
    benchmark::main 
    Threads::Threads
)

# Enforce strict optimization flags for the benchmark target
target_compile_options(bench_monte_carlo PRIVATE 
    -O3 
    -ffast-math 
    -march=native
)

# Register target with CTest
add_test(NAME BenchmarkMonteCarlo COMMAND bench_monte_carlo)