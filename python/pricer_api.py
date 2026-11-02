# python/research/surface_fit.py
import sys
import os
import time

# Ensure the build directory is in the Python path to locate pricing_engine.so
build_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../build/bindings'))
sys.path.insert(0, build_path)

try:
    import pricing_engine as pe
except ImportError as e:
    print(f"Failed to import C++ pricing_engine. Ensure you ran 'make' in the build dir. Error: {e}")
    sys.exit(1)

def run_heston_simulation():
    print("--- High-Performance Heston C++ Pricer ---")
    
    # 1. Define Market & Contract Details
    spot_price = 100.0
    strike_price = 105.0
    risk_free_rate = 0.05
    time_to_maturity = 1.0  # 1 Year
    
    # 2. Instantiate C++ Payoff Object
    call_option = pe.CallPayoff(strike_price)
    
    # 3. Define Heston Parameters (v0, kappa, theta, sigma, rho)
    params = pe.HestonParams(0.04, 2.0, 0.04, 0.3, -0.7)
    
    # 4. Instantiate C++ Pricer
    pricer = pe.HestonPricer(params, risk_free_rate)
    
    # 5. Execute Simulation
    num_paths = 100_000
    num_steps = 252  # Daily steps for 1 year
    
    start_time = time.perf_counter()
    
    # The GIL is released here implicitly if we configure Pybind11 to do so, 
    # executing pure C++ natively.
    price = pricer.price_european_mc(
        spot=spot_price, 
        time_to_maturity=time_to_maturity, 
        num_paths=num_paths, 
        num_time_steps=num_steps, 
        payoff=call_option
    )
    
    end_time = time.perf_counter()
    
    print(f"Option Price:     ${price:.4f}")
    print(f"Simulated Paths:  {num_paths:,}")
    print(f"Time Steps:       {num_steps}")
    print(f"Execution Latency: {(end_time - start_time) * 1000:.2f} ms")

if __name__ == "__main__":
    run_heston_simulation()