#pragma once

namespace hpdf {
namespace constants {

// --- Numeric thresholds ---
constexpr double EPSILON         = 1e-12;
constexpr double SQRT_EPSILON    = 1e-6;
constexpr double TOLERANCE       = 1e-8;
constexpr double PRICE_TOL       = 1e-6;
constexpr double GREEK_TOL       = 1e-5;

// --- Algorithmic limits ---
constexpr int MAX_NEWTON_STEPS   = 32;
constexpr int MAX_BISECT_STEPS   = 64;

// --- Default market parameters ---
constexpr double DEFAULT_R       = 0.01;   // Default risk-free rate
constexpr double DEFAULT_Q       = 0.0;    // Default dividend yield

// --- Compile-time feature toggles ---
constexpr bool ENABLE_FAST_MATH  = true;
constexpr bool ENABLE_SIMD       = false;  // Set true if SIMD code paths are implemented

} // namespace constants
}