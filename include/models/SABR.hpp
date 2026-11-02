// include/models/SABR.hpp
#pragma once

#include <cmath>
#include <stdexcept>

namespace quant::pricing::models {

    struct SABRParams {
        double alpha; // Initial volatility
        double beta;  // CEV elasticity parameter [0, 1]
        double rho;   // Correlation
        double nu;    // Volatility of volatility
    };

    class SABRPricer {
    public:
        explicit SABRPricer(const SABRParams& params) noexcept : m_params(params) {}

        // Hagan's Implied Volatility Expansion (2002)
        [[nodiscard]] double implied_volatility(double F, double K, double T) const {
            if (F <= 0.0 || K <= 0.0 || T <= 0.0) {
                throw std::invalid_argument("Forward, Strike, and Maturity must be positive.");
            }

            const double alpha = m_params.alpha;
            const double beta = m_params.beta;
            const double rho = m_params.rho;
            const double nu = m_params.nu;

            // At-the-money edge case
            if (std::abs(F - K) < 1e-7) {
                const double fkb = std::pow(F, 1.0 - beta);
                const double term1 = alpha / fkb;
                const double term2 = 1.0 + (
                    ((1.0 - beta) * (1.0 - beta) * alpha * alpha) / (24.0 * fkb * fkb) +
                    (rho * beta * nu * alpha) / (4.0 * fkb) +
                    ((2.0 - 3.0 * rho * rho) * nu * nu) / 24.0
                ) * T;
                return term1 * term2;
            }

            const double logFK = std::log(F / K);
            const double fkb = std::pow(F * K, (1.0 - beta) / 2.0);
            
            const double z = (nu / alpha) * fkb * logFK;
            
            // Protect against negative roots in x(z) calculation
            const double root_term = 1.0 - 2.0 * rho * z + z * z;
            if (root_term <= 0.0) return 0.0;
            
            const double xz = std::log((std::sqrt(root_term) + z - rho) / (1.0 - rho));

            const double term1 = alpha / (fkb * (1.0 + ((1.0 - beta) * (1.0 - beta) / 24.0) * logFK * logFK +
                                                 ((1.0 - beta) * (1.0 - beta) * (1.0 - beta) * (1.0 - beta) / 1920.0) * std::pow(logFK, 4)));
            const double term2 = z / xz;
            const double term3 = 1.0 + (
                ((1.0 - beta) * (1.0 - beta) * alpha * alpha) / (24.0 * fkb * fkb) +
                (rho * beta * nu * alpha) / (4.0 * fkb) +
                ((2.0 - 3.0 * rho * rho) * nu * nu) / 24.0
            ) * T;

            return term1 * term2 * term3;
        }

    private:
        SABRParams m_params;
    };

} // namespace quant::pricing::models