# python/tests/test_bindings.py
import sys
import os
import pytest

# Adjust path to find the compiled C++ extension
build_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../build/bindings'))
sys.path.insert(0, build_path)

import pricing_engine as pe

def test_payoff_evaluation():
    """Verify that C++ overridden virtual operators resolve correctly in Python."""
    call = pe.CallPayoff(100.0)
    assert call(110.0) == 10.0, "In-the-money call should return Spot - Strike"
    assert call(90.0) == 0.0, "Out-of-the-money call should return 0.0"
    
    put = pe.PutPayoff(100.0)
    assert put(90.0) == 10.0, "In-the-money put should return Strike - Spot"
    assert put(110.0) == 0.0, "Out-of-the-money put should return 0.0"

def test_heston_struct_binding():
    """Verify that we can read and write to the C++ memory struct natively."""
    params = pe.HestonParams(0.04, 1.5, 0.04, 0.3, -0.5)
    
    assert params.v0 == 0.04
    assert params.kappa == 1.5
    
    # Mutate the struct in Python and verify it holds
    params.rho = -0.8
    assert params.rho == -0.8

def test_heston_mc_execution():
    """Run a small Monte Carlo execution to ensure no segfaults or nan values."""
    params = pe.HestonParams(0.04, 2.0, 0.04, 0.3, -0.7)
    pricer = pe.HestonPricer(params, 0.05)
    call = pe.CallPayoff(100.0)
    
    price = pricer.price_european_mc(100.0, 1.0, 1000, 100, call)
    
    assert price > 0.0, "Option price should be strictly positive."
    assert price < 100.0, "Option price cannot exceed spot price."