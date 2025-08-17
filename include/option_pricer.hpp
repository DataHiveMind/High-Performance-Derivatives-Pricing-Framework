#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include <optional>

namespace hpdf {

// Option type enumeration
enum class OptionType { Call, Put };

// Dividend schedule: vector of (time, amount)
using DividendSchedule = std::vector<std::pair<double, double>>;

// Option specification struct
struct OptionSpec {
    OptionType type;
    double S; // Spot
    double K; // Strike
    double r; // Risk-free rate
    double q; // Dividend yield
    double sigma; // Volatility
    double T; // Time to maturity (years)
    DividendSchedule dividends;

    // Optional: Validate basic invariants
    void validate() const {
        if (S <= 0.0) throw std::invalid_argument("Spot must be positive");
        if (K <= 0.0) throw std::invalid_argument("Strike must be positive");
        if (sigma < 0.0) throw std::invalid_argument("Volatility must be non-negative");
        if (T < 0.0) throw std::invalid_argument("Time to maturity must be non-negative");
        // Additional checks can be added as needed
    }
};

// Greeks container
struct Greeks {
    double delta = 0.0;
    double gamma = 0.0;
    double vega  = 0.0;
    double theta = 0.0;
    double rho   = 0.0;
    // Add more if needed
};

// Price result container
struct PriceResult {
    double price = 0.0;
    Greeks greeks;
};

// Abstract base class for option pricers
struct OptionPricer {
    virtual ~OptionPricer() = default;

    // Compute price and greeks
    virtual PriceResult price(const OptionSpec& spec) const = 0;

    // Compute only greeks (optional, can default to price().greeks)
    virtual Greeks greeks(const OptionSpec& spec) const {
        return price(spec).greeks;
    }

    // Validate input spec (optional override)
    virtual void validate(const OptionSpec& spec) const {
        spec.validate();
    }
};