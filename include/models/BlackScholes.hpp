// include/models/BlackScholes.hpp
#pragma once

#include "core/payoff.hpp"
#include <cmath>
#include <numbers>

namespace quant::pricing::models {

    class BlackScholes {
    public:
        // Cumulative Distribution Function for Standard Normal
        [[nodiscard]] static double norm_cdf(double x) noexcept {
            return 0.5 * std::erfc(-x * std::numbers::inv_sqrt2);
        }

        // Analytical European Option Pricing
        [[nodiscard]] static double price_european(
            double S, double K, double r, double v, double T, bool is_call) noexcept 
        {
            if (T <= 0.0) {
                return is_call ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
            }

            const double d1 = (std::log(S / K) + (r + 0.5 * v * v) * T) / (v * std::sqrt(T));
            const double d2 = d1 - v * std::sqrt(T);

            if (is_call) {
                return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
            } else {
                return K * std::exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1);
            }
        }
    };

} // namespace quant::pricing::models