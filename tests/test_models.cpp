// tests/test_models.cpp
#include <gtest/gtest.h>
#include "models/Heston.hpp"
#include "models/SABR.hpp"
#include "core/payoff.hpp"

using namespace quant::pricing;

TEST(HestonModelTest, MonteCarloConvergence) {
    // Parameters: v0, kappa, theta, sigma, rho
    models::HestonParams params{0.04, 2.0, 0.04, 0.1, 0.0}; 
    models::HestonPricer pricer(params, 0.05);
    core::CallPayoff call(100.0);

    // If volatility of volatility (sigma) is low, and v0 == theta, 
    // Heston price should closely approximate Black-Scholes with v=0.20
    double mc_price = pricer.price_european_mc(100.0, 1.0, 50000, 100, call);
    
    // BS Price for S=100, K=100, r=0.05, v=0.20, T=1 is ~10.45
    EXPECT_NEAR(mc_price, 10.45, 0.20); 
    EXPECT_GT(mc_price, 0.0); // Strict non-negativity bound
}

TEST(SABRModelTest, HaganImpliedVolatilityATM) {
    // Parameters: alpha, beta, rho, nu
    models::SABRParams params{0.20, 1.0, -0.5, 0.4};
    models::SABRPricer pricer(params);
    
    // Forward = Strike = 100
    double implied_vol = pricer.implied_volatility(100.0, 100.0, 1.0);
    
    // ATM vol should be strictly positive and near the initial alpha parameter
    EXPECT_GT(implied_vol, 0.0);
    EXPECT_NEAR(implied_vol, 0.20, 0.05);
}

TEST(SABRModelTest, InvalidInputsThrowOrZero) {
    models::SABRParams params{0.20, 1.0, -0.5, 0.4};
    models::SABRPricer pricer(params);
    
    EXPECT_THROW(pricer.implied_volatility(-100.0, 100.0, 1.0), std::invalid_argument);
    EXPECT_THROW(pricer.implied_volatility(100.0, 100.0, -1.0), std::invalid_argument);
}