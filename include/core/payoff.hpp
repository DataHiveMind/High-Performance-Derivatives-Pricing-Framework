// include/core/payoff.hpp
#pragma once

#include <algorithm>
#include <concepts>

namespace quant::pricing::core {

    // Abstract Base Class for Payoffs
    class IPayoff {
    public:
        virtual ~IPayoff() = default;
        
        // Pure virtual function overloaded for functor evaluation
        [[nodiscard]] virtual double operator()(double spot) const noexcept = 0;
    };

    // Standard European Call Option
    class CallPayoff final : public IPayoff {
    public:
        explicit CallPayoff(double strike) noexcept : m_strike(strike) {}

        [[nodiscard]] double operator()(double spot) const noexcept override {
            return std::max(spot - m_strike, 0.0);
        }

    private:
        double m_strike;
    };

    // Standard European Put Option
    class PutPayoff final : public IPayoff {
    public:
        explicit PutPayoff(double strike) noexcept : m_strike(strike) {}

        [[nodiscard]] double operator()(double spot) const noexcept override {
            return std::max(m_strike - spot, 0.0);
        }

    private:
        double m_strike;
    };

} // namespace quant::pricing::core