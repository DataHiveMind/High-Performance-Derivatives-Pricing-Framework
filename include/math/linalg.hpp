// include/math/linalg.hpp
#pragma once

#include <vector>
#include <cmath>
#include <array>
#include <stdexcept>

namespace quant::pricing::math {

    // SIMD-aligned memory allocator (AVX-512 requires 64-byte alignment)
    template <typename T, std::size_t Alignment = 64>
    struct AlignedAllocator {
        using value_type = T;

        AlignedAllocator() noexcept = default;
        template <class U> AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

        T* allocate(std::size_t n) {
            void* ptr = nullptr;
            if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0) {
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        void deallocate(T* ptr, std::size_t) noexcept {
            free(ptr);
        }
    };

    template <typename T>
    using aligned_vector = std::vector<T, AlignedAllocator<T>>;

    // 2x2 Cholesky Decomposition for correlating random variables inline
    [[nodiscard]] constexpr std::array<std::array<double, 2>, 2> 
    cholesky_2x2(double rho) noexcept {
        return {{
            {1.0, 0.0},
            {rho, std::sqrt(1.0 - rho * rho)}
        }};
    }

} // namespace quant::pricing::math