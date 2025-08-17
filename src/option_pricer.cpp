#include "option_pricer.hpp"
#include <stdexcept>
#include <sstream>

namespace hpdf {

// Common validation utility (can be extended for logging or diagnostics)
void validate_option_spec(const OptionSpec& spec) {
    try {
        spec.validate();
    } catch (const std::exception& ex) {
        std::ostringstream oss;
        oss << "OptionSpec validation failed: " << ex.what();
        throw std::invalid_argument(oss.str());
    }
}

// Adapter: Validate and price using any OptionPricer
PriceResult price_with_validation(const OptionPricer& pricer, const OptionSpec& spec) {
    validate_option_spec(spec);
    return pricer.price(spec);
}

} // namespace