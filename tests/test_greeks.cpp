// tests/test_greeks.cpp
#include <gtest/gtest.h>
#include "models/BlackScholes.hpp"
#include <cmath>

using namespace quant::pricing::models;

class GreeksTest : public ::testing::Test {
protected:
    const double S = 100.0;
    const double K = 100.0;
    const double r = 0.05;
    const double v = 0.20;
    const double T = 1.0;
    const double dS = 0.01; // 1 cent bump
    const double dv = 0.001; // 10 bps vol bump
};

TEST_F(GreeksTest, CallDeltaCentralDifference) {
    // Bump up and down
    double price_up = BlackScholes::price_european(S + dS, K, r, v, T, true);
    double price_down = BlackScholes::price_european(S - dS, K, r, v, T, true);
    
    double numerical_delta = (price_up - price_down) / (2.0 * dS);
    
    // Analytical Delta for ATM Call is roughly 0.5 to 0.6
    EXPECT_NEAR(numerical_delta, 0.6368, 1e-3);
}

TEST_F(GreeksTest, CallGammaCentralDifference) {
    double price_up = BlackScholes::price_european(S + dS, K, r, v, T, true);
    double price_mid = BlackScholes::price_european(S, K, r, v, T, true);
    double price_down = BlackScholes::price_european(S - dS, K, r, v, T, true);
    
    double numerical_gamma = (price_up - 2.0 * price_mid + price_down) / (dS * dS);
    
    EXPECT_NEAR(numerical_gamma, 0.0187, 1e-3);
}

TEST_F(GreeksTest, CallVegaForwardDifference) {
    double price_up = BlackScholes::price_european(S, K, r, v + dv, T, true);
    double price_mid = BlackScholes::price_european(S, K, r, v, T, true);
    
    double numerical_vega = (price_up - price_mid) / dv;
    
    // Standard Vega is quoted per 1% change, but we calculated raw partial derivative
    EXPECT_NEAR(numerical_vega, 37.524, 1e-2);
}