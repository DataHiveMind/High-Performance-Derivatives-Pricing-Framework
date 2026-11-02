// include/models/Heston.hpp
#pragma once

#include "core/payoff.hpp"
#include <vector>
#include <memory>
#include <random>

namespace quant::pricing::models {

    // Heston Model Parameters Structure
    struct HestonParams {
        double v0;    // Initial variance
        double kappa; // Mean reversion speed
        double theta; // Long-term mean variance
        double sigma; // Volatility of volatility
        double rho;   // Correlation between asset price and volatility Wiener processes
    };

    class HestonPricer {
    public:
        HestonParams params;
        double risk_free_rate;

        explicit HestonPricer(const HestonParams& p, double r) noexcept 
            : params(p), risk_free_rate(r) {}

        // Monte Carlo pricing engine utilizing Euler-Maruyama discretization
        [[nodiscard]] double price_european_mc(
            double spot, 
            double time_to_maturity, 
            size_t num_paths, 
            size_t num_time_steps,
            const core::IPayoff& payoff
        ) const;
    };

} // namespace quant::pricing::models