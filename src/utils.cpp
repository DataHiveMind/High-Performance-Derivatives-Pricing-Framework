#include "utils.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <mutex>

namespace hpdf {
namespace utils {

// --- Math: Normal PDF/CDF (fast approximations) ---

// Standard normal PDF (thread-safe, stateless)
double norm_pdf(double x) {
    static constexpr double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

// Standard normal CDF (erf-based, accurate to ~1e-7)
double norm_cdf(double x) {
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

// --- I/O: Minimal CSV loader ---

std::vector<std::vector<std::string>> load_csv(const std::string& path, char sep) {
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

// --- Timing: Scoped timer, now_ms ---

uint64_t now_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

// --- Random: Thread-safe PRNG seeding ---

std::mt19937 seeded_prng() {
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    return std::mt19937(seq);
}

} // namespace utils
} //