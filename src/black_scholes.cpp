#include "black_scholes.hpp"
#include <cmath>
#include <algorithm>

namespace hpdf {
namespace black_scholes {

// --- Internal: Stable d1/d2 computation (implementation) ---
void compute_d1_d2(double S, double K, double r, double q, double sigma, double T, double& d1, double& d2) {
    if (sigma <= 0.0 || T <= 0.0) {
        d1 = d2 = 0.0;
        return;
    }
    double vsqrt = sigma * std::sqrt(T);
    // Use log1p/exp_m1 for better accuracy if S ≈ K
    d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / vsqrt;
    d2 = d1 - vsqrt;
}

// --- Batch pricing (vectorized loop for batch mode) ---
void price_european_batch(
    const std::vector<OptionSpec>& specs,
    std::vector<PriceResult>& results)
{
    results.resize(specs.size());
    for (size_t i = 0; i < specs.size(); ++i) {
        const auto& spec = specs[i];
        if (spec.type == OptionType::Call) {
            results[i].price = price_european_call(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
            results[i].greeks = greeks_call(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
        } else {
            results[i].price = price_european_put(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
            results[i].greeks = greeks_put(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
        }
    }
}

// --- Theta conventions ---
// By default, theta is per calendar day. For trading day, scale by ~1/252.

double theta_trading_day(const Greeks& g) {
    return g.theta / 252.0;
}

double theta_calendar_day(const Greeks& g) {
    return g.theta / 365.0;
}

} // namespace black_scholes
} // namespace hpdf