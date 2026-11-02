// include/math/rng.hpp
#pragma once

#include <cstdint>
#include <limits>
#include <cmath>
#include <numbers>

namespace quant::pricing::math {

    // Xoshiro256** PRNG for extreme throughput
    class Xoshiro256StarStar {
    public:
        explicit Xoshiro256StarStar(uint64_t seed) noexcept {
            // SplitMix64 initialization
            uint64_t z = seed;
            for (int i = 0; i < 4; ++i) {
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
                z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
                s[i] = z ^ (z >> 31);
            }
        }

        [[nodiscard]] uint64_t operator()() noexcept {
            const uint64_t result = rotl(s[1] * 5, 7) * 9;
            const uint64_t t = s[1] << 17;
            s[2] ^= s[0];
            s[3] ^= s[1];
            s[1] ^= s[2];
            s[0] ^= s[3];
            s[2] ^= t;
            s[3] = rotl(s[3], 45);
            return result;
        }

        // Generate uniform double in [0, 1)
        [[nodiscard]] double next_double() noexcept {
            return (operator()() >> 11) * 0x1.0p-53;
        }

        // Generate Standard Normal N(0,1) via Box-Muller
        [[nodiscard]] double next_normal() noexcept {
            double u1 = next_double();
            double u2 = next_double();
            // Protect against log(0)
            u1 = (u1 < 1e-15) ? 1e-15 : u1;
            return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * std::numbers::pi * u2);
        }

    private:
        uint64_t s[4];
        static inline uint64_t rotl(const uint64_t x, int k) noexcept {
            return (x << k) | (x >> (64 - k));
        }
    };

} // namespace quant::pricing::math