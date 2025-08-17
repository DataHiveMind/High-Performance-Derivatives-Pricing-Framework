#pragma once
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>
#include <functional>
#include <iostream>

namespace hpdf {
namespace utils {

// --- Math: Normal PDF/CDF, clamp, stable exp/log ---

inline double norm_pdf(double x) {
    static constexpr double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

inline double norm_cdf(double x) {
    // Using std::erf for numerical stability
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

template<typename T>
constexpr T clamp(T v, T lo, T hi) {
    return std::max(lo, std::min(v, hi));
}

// Stable exp/log helpers
inline double safe_exp(double x) {
    constexpr double max_exp = 700.0; // ~exp(700) is max double
    constexpr double min_exp = -700.0;
    return std::exp(clamp(x, min_exp, max_exp));
}

inline double safe_log(double x) {
    constexpr double min_x = 1e-300;
    return std::log(std::max(x, min_x));
}

// --- I/O: Minimal CSV loader, parameter parsing ---

// Loads a CSV file into a vector of vector<string>
inline std::vector<std::vector<std::string>> load_csv(const std::string& path, char sep = ',') {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open CSV: " + path);
    std::vector<std::vector<std::string>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::istringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, sep)) {
            row.push_back(cell);
        }
        data.push_back(row);
    }
    return data;
}

// Parse parameter from string with fallback
template<typename T>
T parse_param(const std::string& s, T fallback) {
    std::istringstream iss(s);
    T val;
    if (!(iss >> val)) return fallback;
    return val;
}

// --- Timing: Scoped timer, steady_clock wrappers ---

class ScopedTimer {
    std::string label_;
    std::chrono::steady_clock::time_point start_;
public:
    explicit ScopedTimer(const std::string& label)
        : label_(label), start_(std::chrono::steady_clock::now()) {}
    ~ScopedTimer() {
        auto end = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
        std::cerr << "[Timer] " << label_ << ": " << ms << " ms\n";
    }
};

inline uint64_t now_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

// --- Random: PRNG seeding ---

inline std::mt19937 seeded_prng() {
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    return std::mt19937(seq);
}

} // namespace utils
} //