// src/pde_solvers.cpp
#include <vector>
#include <stdexcept>

namespace quant::pricing::pde {

    // Solves the tridiagonal system matrix equation Ax = d
    // a: sub-diagonal, b: main diagonal, c: super-diagonal, d: right-hand side
    void solve_tridiagonal_system(
        const std::vector<double>& a, 
        std::vector<double>& b, 
        const std::vector<double>& c, 
        std::vector<double>& d, 
        std::vector<double>& x) 
    {
        const size_t n = d.size();
        if (n == 0) return;
        if (a.size() != n || b.size() != n || c.size() != n || x.size() != n) {
            throw std::invalid_argument("Tridiagonal solver vectors must be of equal size.");
        }

        // Forward elimination phase
        for (size_t i = 1; i < n; ++i) {
            const double m = a[i] / b[i - 1];
            b[i] -= m * c[i - 1];
            d[i] -= m * d[i - 1];
        }

        // Backward substitution phase
        x[n - 1] = d[n - 1] / b[n - 1];
        for (size_t i = n - 1; i-- > 0;) {
            x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
        }
    }

} // namespace quant::pricing::pde