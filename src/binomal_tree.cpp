#include "binomial_tree.hpp"
#include <cmath>
#include <algorithm>

namespace hpdf {

// --- Internal: Factor generation and risk-neutral probability ---
void compute_binomial_factors(
    double r, double q, double sigma, double dt,
    BinomialScheme scheme,
    double& u, double& d, double& p)
{
    if (scheme == BinomialScheme::CRR) {
        u = std::exp(sigma * std::sqrt(dt));
        d = 1.0 / u;
        p = (std::exp((r - q) * dt) - d) / (u - d);
    } else { // Jarrow–Rudd
        u = std::exp((r - q - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt));
        d = std::exp((r - q - 0.5 * sigma * sigma) * dt - sigma * std::sqrt(dt));
        p = 0.5;
    }
    p = utils::clamp(p, 0.0, 1.0);
}

// --- Core lattice pricer (single buffer, branch-minimized) ---
PriceResult price_binomial_impl(
    const OptionSpec& spec,
    int steps,
    BinomialScheme scheme)
{
    const double S = spec.S;
    const double K = spec.K;
    const double r = spec.r;
    const double q = spec.q;
    const double sigma = spec.sigma;
    const double T = spec.T;
    const bool is_call = (spec.type == OptionType::Call);
    const bool is_american = !spec.dividends.empty(); // crude: treat non-empty dividends as American

    const double dt = T / steps;
    double u, d, p;
    compute_binomial_factors(r, q, sigma, dt, scheme, u, d, p);

    const double disc = std::exp(-r * dt);

    std::vector<double> values(steps + 1);

    // Terminal payoffs
    for (int i = 0; i <= steps; ++i) {
        double ST = S * std::pow(u, steps - i) * std::pow(d, i);
        values[i] = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
    }

    // Backward induction
    for (int step = steps - 1; step >= 0; --step) {
        #pragma omp parallel for if(steps > 64)
        for (int i = 0; i <= step; ++i) {
            double cont = disc * (p * values[i] + (1.0 - p) * values[i + 1]);
            if (is_american) {
                double ST = S * std::pow(u, step - i) * std::pow(d, i);
                double ex = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
                values[i] = std::max(cont, ex);
            } else {
                values[i] = cont;
            }
        }
    }

    PriceResult result;
    result.price = values[0];
    // Greeks via bump-and-reprice (finite difference)
    // For performance, user should call BinomialTreePricer::greeks if needed.
    return result;
}

}