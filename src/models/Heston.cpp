// src/models/Heston.cpp
#include "models/Heston.hpp"
#include "math/rng.hpp"
#include "math/linalg.hpp"
#include <cmath>
#include <algorithm>

namespace quant::pricing::models {

    double HestonPricer::price_european_mc(
        double spot, 
        double time_to_maturity, 
        size_t num_paths, 
        size_t num_time_steps,
        const core::IPayoff& payoff
    ) const {
        if (time_to_maturity <= 0.0) {
            return payoff(spot);
        }

        const double dt = time_to_maturity / static_cast<double>(num_time_steps);
        const double sqrt_dt = std::sqrt(dt);
        
        // Pre-compute Cholesky decomposition for correlated Brownian motions
        const auto cholesky = math::cholesky_2x2(params.rho);
        const double rho_c11 = cholesky[0][0]; // 1.0
        const double rho_c21 = cholesky[1][0]; // rho
        const double rho_c22 = cholesky[1][1]; // sqrt(1 - rho^2)

        // Initialize high-throughput RNG (Seeded statically for deterministic benchmarking)
        math::Xoshiro256StarStar rng(42ULL);

        double sum_payoffs = 0.0;

        // Path generation loop
        for (size_t i = 0; i < num_paths; ++i) {
            double s_t = spot;
            double v_t = params.v0;

            for (size_t t = 0; t < num_time_steps; ++t) {
                // Generate independent standard normals
                const double z1 = rng.next_normal();
                const double z2 = rng.next_normal();

                // Correlate the random variables
                const double dw1 = z1 * sqrt_dt;
                const double dw2 = (rho_c21 * z1 + rho_c22 * z2) * sqrt_dt;

                // Ensure variance remains non-negative (Full Truncation scheme)
                const double v_t_plus = std::max(v_t, 0.0);
                const double sqrt_v_t = std::sqrt(v_t_plus);

                // Euler-Maruyama step for Asset Price (using log-process to prevent negative prices)
                s_t *= std::exp((risk_free_rate - 0.5 * v_t_plus) * dt + sqrt_v_t * dw1);

                // Euler-Maruyama step for Variance Process
                v_t += params.kappa * (params.theta - v_t_plus) * dt + params.sigma * sqrt_v_t * dw2;
            }

            sum_payoffs += payoff(s_t);
        }

        // Discount the expected payoff to present value
        return (sum_payoffs / static_cast<double>(num_paths)) * std::exp(-risk_free_rate * time_to_maturity);
    }

} // namespace quant::pricing::models