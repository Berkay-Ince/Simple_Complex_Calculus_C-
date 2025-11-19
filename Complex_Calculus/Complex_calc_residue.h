#pragma once
#include "Complex.h"
#include "Complex_calc.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <numbers>
#include <functional> 

namespace cc {

    // ============================================
    // 1) Where function is not analytic in the given mesh
    // ============================================

    template<typename F>
    std::vector<Complex>
    finding_singularity(F f,
                        double x_min, double x_max,
                        double y_min, double y_max,
                        int N = 2000)
    {
        double delta_x = std::abs(x_max - x_min) / N;
        double delta_y = std::abs(y_max - y_min) / N;

        std::vector<Complex> singularity_set;
        singularity_set.reserve(N);  // at most N points

        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
                double x = x_min + kx * delta_x;
                double y = y_min + ky * delta_y;
                Complex z{x, y};

                if (!cc::is_analytic(f, z)) {   // you already have this
                    singularity_set.push_back(z);
                }
            }
        }
        return singularity_set;
    }

    // ============================================
    // 2) Finding Degree of Poles (for each singularity)
    // ============================================

    template<typename F>
    std::vector<int>
    finding_degree_of_poles(F f,
                            double x_min, double x_max,
                            double y_min, double y_max,
                            int N = 2000,
                            int max_degree = 10)
    {
        std::vector<Complex> set =
            finding_singularity(f, x_min, x_max, y_min, y_max, N);

        std::vector<int> degrees;
        degrees.reserve(set.size());

        for (std::size_t k = 0; k < set.size(); ++k) {
            Complex z0 = set[k];
            int d = -1;   // -1 = "could not determine"

            for (int m = 1; m <= max_degree; ++m) {

                // g_m(z) = (z - z0)^m * f(z)
                auto g_m = [=](Complex z) -> Complex {
                    return cc::pow_int(z - z0, m) * f(z);
                }; // lamda functions are a bit confusing they seems to me it is a method while defining functions inside the nested structures, but it has also some language source reasons

                if (cc::is_analytic(g_m, z0)) {
                    d = m;
                    break;
                }
            }

            degrees.push_back(d);
        }

        return degrees;
    }

    // ============================================
    // 3) Residue of a Function at a given point z0, with known order m
    // ============================================
    template<typename F>
    Complex residue_at(F f, Complex z0, int m, double h = 1e-5)
    {
        if (m <= 0) {
            return {0.0, 0.0};  // not a pole
        }

        // g_m(z) = (z - z0)^m * f(z)
        auto g_m = [=](Complex z) -> Complex {
            return cc::pow_int(z - z0, m) * f(z);  // be sure pow_int takes (Complex,int)
        };

        // Simple pole: m = 1
        // g_1(z) = (z - z0) f(z) = a_{-1} + ... ⇒ g_1(z0) = residue
        if (m == 1) {
            Complex near_z0{ z0.x + h, z0.y };
            return g_m(near_z0);  // ≈ residue (numerical)
        }

        // m > 1: Res = g_m^{(m-1)}(z0) / (m-1)!
        std::function<Complex(Complex)> current = g_m;

        for (int r = 0; r < m - 1; ++r) {
            auto previous = current;
            current = [=](Complex z) -> Complex {
                return cc::complex_derivative(previous, z, h);
            };
        }

        // current ≈ g_m^{(m-1)}
        Complex g_m_deriv = current(z0);
        double fact = std::tgamma(static_cast<double>(m)); // (m-1)! = Γ(m)
        return g_m_deriv / fact;
    }

    // ============================================
    // 4) Residues of f on a mesh (vector of residues)
    // ============================================

    template<typename F>
    std::vector<Complex>
    residues(F f,double x_min, double x_max, double y_min, double y_max,int meshN = 2000,int max_degree = 10)
    {
        std::vector<Complex> set_1 =
            finding_singularity(f, x_min, x_max, y_min, y_max, meshN);

        std::vector<int> set_2 =
            finding_degree_of_poles(f, x_min, x_max, y_min, y_max, meshN, max_degree);

        std::vector<Complex> result;
        result.reserve(set_1.size());

        for (std::size_t k = 0; k < set_1.size(); ++k) {
            Complex z0 = set_1[k];
            int m      = set_2[k];

            Complex res = residue_at(f, z0, m);
            result.push_back(res);
        }

        return result;
    }

    // ============================================
    // 5) Residue Integral (contour integral) via residue theorem
    // ============================================

    template<typename F>
    Complex residue_integral(F f,
                             double x_min, double x_max,
                             double y_min, double y_max,
                             int meshN      = 2000,
                             int max_degree = 10)
    {
        std::vector<Complex> set =
            cc::residues(f, x_min, x_max, y_min, y_max, meshN, max_degree);

        Complex final{0.0, 0.0};
        for (std::size_t k = 0; k < set.size(); ++k) {
            final = final + set[k];  // or final += set[k] if you have operator+=
        }

        // 2πi * sum residues
        Complex coeff{0.0, 2.0};
        return coeff * final;
    }

} // namespace cc
