// include/core/instruments.hpp
#pragma once

#include "core/payoff.hpp"
#include <memory>
#include <stdexcept>

namespace quant::pricing::core {

    enum class OptionStyle { European, American, Bermudan, Asian };

    class Instrument {
    public:
        virtual ~Instrument() = default;
        
        [[nodiscard]] virtual OptionStyle get_style() const noexcept = 0;
        [[nodiscard]] virtual double get_time_to_maturity() const noexcept = 0;
        [[nodiscard]] virtual const IPayoff& get_payoff() const noexcept = 0;
    };

    class EuropeanOption final : public Instrument {
    public:
        EuropeanOption(std::shared_ptr<IPayoff> payoff, double time_to_maturity)
            : m_payoff(std::move(payoff)), m_ttm(time_to_maturity) {
            if (m_ttm < 0.0) throw std::invalid_argument("Maturity cannot be negative.");
            if (!m_payoff) throw std::invalid_argument("Payoff cannot be null.");
        }

        [[nodiscard]] OptionStyle get_style() const noexcept override { return OptionStyle::European; }
        [[nodiscard]] double get_time_to_maturity() const noexcept override { return m_ttm; }
        [[nodiscard]] const IPayoff& get_payoff() const noexcept override { return *m_payoff; }

    private:
        std::shared_ptr<IPayoff> m_payoff;
        double m_ttm;
    };

} // namespace quant::pricing::core