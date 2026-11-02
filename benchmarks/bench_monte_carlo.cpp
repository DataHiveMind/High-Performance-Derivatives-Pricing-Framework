// benchmarks/bench_monte_carlo.cpp
#include <benchmark/benchmark.h>
#include <vector>
#include <random>

// Mock inclusion of your core engine
// #include "models/Heston.hpp"
// #include "math/rng.hpp"

// Placeholder for the actual vectorized Monte Carlo function
static double price_european_mc(double S, double K, double r, double v, double T, size_t num_paths) {
    // Simulated math overhead
    double sum_payoffs = 0.0;
    for (size_t i = 0; i < num_paths; ++i) {
        // SIMD generation would replace this in the real framework
        double ST = S * std::exp((r - 0.5 * v * v) * T + v * std::sqrt(T) * 0.1); 
        sum_payoffs += std::max(ST - K, 0.0);
    }
    return (sum_payoffs / num_paths) * std::exp(-r * T);
}

// Benchmark: European Call Monte Carlo Engine
static void BM_MonteCarlo_EuropeanCall(benchmark::State& state) {
    // Option parameters
    const double S = 100.0;
    const double K = 100.0;
    const double r = 0.05;
    const double v = 0.20;
    const double T = 1.0;
    
    // Number of paths is passed dynamically via state.range(0)
    const size_t num_paths = state.range(0);

    for (auto _ : state) {
        // DoNotOptimize ensures the compiler doesn't delete the math computation
        double price = price_european_mc(S, K, r, v, T, num_paths);
        benchmark::DoNotOptimize(price);
        
        // ClobberMemory forces memory writes to complete, preventing cache-tricking
        benchmark::ClobberMemory();
    }
    
    // Metrics reporting: Allows you to see "Paths per Second" in the output
    state.SetItemsProcessed(state.iterations() * num_paths);
}

// Register the function as a benchmark
// Test across different magnitudes of paths: 10k, 100k, 1 Million
BENCHMARK(BM_MonteCarlo_EuropeanCall)->RangeMultiplier(10)->Range(10000, 1000000);

BENCHMARK_MAIN();