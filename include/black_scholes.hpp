#pragma once
#include "option_pricer.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include <cmath>
#include <algorithm>

namespace hpdf {

// --- Analytic Black-Scholes(-Merton) pricer for European options ---
// Supports continuous dividend yield q. For discrete dividends, see notes below.

namespace black_scholes {

// --- Internal: Stable d1/d2 computation ---
inline void compute_d1_d2(double S, double K, double r, double q, double sigma, double T, double& d1, double& d2) {
    if (sigma <= 0.0 || T <= 0.0) {
        d1 = d2 = 0.0;
        return;
    }
    double vsqrt = sigma * std::sqrt(T);
    d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / vsqrt;
    d2 = d1 - vsqrt;
}

// --- Price formulas (Merton: with continuous dividend yield) ---
inline double price_european_call(double S, double K, double r, double q, double sigma, double T) {
    double d1, d2;
    compute_d1_d2(S, K, r, q, sigma, T, d1, d2);
    double df_r = std::exp(-r * T);
    double df_q = std::exp(-q * T);
    return S * df_q * utils::norm_cdf(d1) - K * df_r * utils::norm_cdf(d2);
}

inline double price_european_put(double S, double K, double r, double q, double sigma, double T) {
    double d1, d2;
    compute_d1_d2(S, K, r, q, sigma, T, d1, d2);
    double df_r = std::exp(-r * T);
    double df_q = std::exp(-q * T);
    return K * df_r * utils::norm_cdf(-d2) - S * df_q * utils::norm_cdf(-d1);
}

// --- Greeks (analytic) ---
inline Greeks greeks_call(double S, double K, double r, double q, double sigma, double T) {
    Greeks g;
    double d1, d2;
    compute_d1_d2(S, K, r, q, sigma, T, d1, d2);
    double df_r = std::exp(-r * T);
    double df_q = std::exp(-q * T);
    double n_d1 = utils::norm_pdf(d1);

    g.delta = df_q * utils::norm_cdf(d1);
    g.gamma = df_q * n_d1 / (S * sigma * std::sqrt(T));
    g.vega  = S * df_q * n_d1 * std::sqrt(T);
    g.theta = - (S * sigma * df_q * n_d1) / (2 * std::sqrt(T))
              - r * K * df_r * utils::norm_cdf(d2)
              + q * S * df_q * utils::norm_cdf(d1);
    g.rho   = K * T * df_r * utils::norm_cdf(d2);
    return g;
}

inline Greeks greeks_put(double S, double K, double r, double q, double sigma, double T) {
    Greeks g;
    double d1, d2;
    compute_d1_d2(S, K, r, q, sigma, T, d1, d2);
    double df_r = std::exp(-r * T);
    double df_q = std::exp(-q * T);
    double n_d1 = utils::norm_pdf(d1);

    g.delta = -df_q * utils::norm_cdf(-d1);
    g.gamma = df_q * n_d1 / (S * sigma * std::sqrt(T));
    g.vega  = S * df_q * n_d1 * std::sqrt(T);
    g.theta = - (S * sigma * df_q * n_d1) / (2 * std::sqrt(T))
              + r * K * df_r * utils::norm_cdf(-d2)
              - q * S * df_q * utils::norm_cdf(-d1);
    g.rho   = -K * T * df_r * utils::norm_cdf(-d2);
    return g;
}

// --- OptionPricer implementation ---
struct BlackScholesPricer : public OptionPricer {
    PriceResult price(const OptionSpec& spec) const override {
        validate(spec);
        PriceResult res;
        if (spec.type == OptionType::Call) {
            res.price = price_european_call(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
            res.greeks = greeks_call(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
        } else {
            res.price = price_european_put(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
            res.greeks = greeks_put(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
        }
        return res;
    }
    Greeks greeks(const OptionSpec& spec) const override {
        if (spec.type == OptionType::Call)
            return greeks_call(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
        else
            return greeks_put(spec.S, spec.K, spec.r, spec.q, spec.sigma, spec.T);
    }
    void validate(const OptionSpec& spec) const override {
        spec.validate();
        if (!spec.dividends.empty())
            throw std::invalid_argument("Black-Scholes: Discrete dividends not supported (use q for continuous yield)");
    }
};

// --- Notes ---
// For discrete dividends, analytic pricing is not closed-form. Approximations include:
//  - Adjusting S for PV(dividends) and using q=0
//  - Using binomial trees or finite difference methods

} // namespace black_scholes
} //