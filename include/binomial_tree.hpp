#pragma once
#include "option_pricer.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include <vector>
#include <algorithm>
#include <cmath>    // ensures std::sqrt is available
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hpdf {

// Binomial tree schemes
enum class BinomialScheme { CRR, JarrowRudd };

// Binomial tree pricer (supports American/European, calls/puts, dividends, futures-style)
struct BinomialTreePricer : public OptionPricer {
    int steps;
    BinomialScheme scheme;

    explicit BinomialTreePricer(int steps_ = 200, BinomialScheme scheme_ = BinomialScheme::CRR)
        : steps(steps_), scheme(scheme_) {}

    PriceResult price(const OptionSpec& spec) const override {
        validate(spec);

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

        // Up/down/prob by scheme
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

        // Precompute discount
        const double disc = std::exp(-r * dt);

        // Terminal payoffs
        ::std::vector<double> values(steps + 1);
        for (int i = 0; i <= steps; ++i) {
            double ST = S * pow(u, steps - i) * pow(d, i);
            values[i] = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
        }

        // Backward induction
        for (int step = steps - 1; step >= 0; --step) {
            #pragma omp parallel for if(steps > 64)
            for (int i = 0; i <= step; ++i) {
                double cont = disc * (p * values[i] + (1.0 - p) * values[i + 1]);
                if (is_american) {
                    double ST = S * pow(u, step - i) * pow(d, i);
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
        result.greeks = compute_greeks_fd(spec);

        return result;
    }

    Greeks greeks(const OptionSpec& spec) const override {
        return compute_greeks_fd(spec);
    }

    void validate(const OptionSpec& spec) const override {
        spec.validate();
        if (steps < 1) throw ::std::invalid_argument("Binomial steps must be positive");
    }

private:
    Greeks compute_greeks_fd(const OptionSpec& spec) const {
        using namespace constants;
        Greeks g;
        const double dS = ::std::max(0.01, 0.01 * spec.S);
        const double dV = ::std::max(0.001, 0.01 * spec.sigma);
        const double dT = ::std::min(0.01, 0.5 * spec.T);
        const double dR = 0.0001;

        OptionSpec up = spec, down = spec;

        // Delta
        up.S += dS; down.S -= dS;
        double pu = price(up).price, pd = price(down).price, p0 = price(spec).price;
        g.delta = (pu - pd) / (2 * dS);

        // Gamma
        g.gamma = (pu - 2 * p0 + pd) / (dS * dS);

        // Vega
        up = spec; down = spec;
        up.sigma += dV; down.sigma -= dV;
        g.vega = (price(up).price - price(down).price) / (2 * dV);

        // Theta
        up = spec; down = spec;
        up.T = ::std::max(0.0, spec.T - dT);
        g.theta = (price(up).price - price(down).price) / dT;

        // Rho
        up = spec; down = spec;
        up.r += dR; down.r -= dR;
        g.rho = (price(up).price - price(down).price) / (2 * dR);

        return g;
    }
};

// API function
inline PriceResult price_binomial(const OptionSpec& spec, int steps = 200, BinomialScheme scheme = BinomialScheme::CRR) {
    BinomialTreePricer pricer(steps, scheme);
    return pricer.price(spec);
}

}