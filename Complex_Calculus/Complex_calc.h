#pragma once
#include "Complex.h"
#include <stdexcept>
#include "cmath"
using std::runtime_error;

namespace cc {
    
    // Analyticity of Function

    template<typename F>
    bool is_analytic(F f, Complex z0, double h= 1e-4, double tol = 1e-3){
        
        double x = z0.x;
        double y = z0.y;

        double ux = (f({x+h,y}).x-f({x-h,y}).x)/(2*h); // f({}).x = u
        double uy = (f({x,y+h}).x-f({x,y-h}).x)/(2*h);
        double vx = (f({x+h,y}).y-f({x-h,y}).y)/(2*h); // f({}).y = v
        double vy = (f({x,y+h}).y-f({x,y-h}).y)/(2*h);

        bool cr1 = std::abs(ux-vy) < tol;
        bool cr2 = std::abs(uy+vx) < tol;

        return cr1 && cr2;
    }
    
    // Complex Derivative
    
    template<typename F>
    Complex complex_derivative(F f,Complex z0, double h = 1e-6){
        
        if (!is_analytic(f,z0)){
            throw runtime_error("Function is not analytic!");
        }

        Complex hC{h,0};
        return (f(z0+hC)-f(z0-hC)) / Complex{2*h,0};   
    }

    // Parametrized Curve

    Complex complex_curve(Complex z0, Complex z1, double t){
        double gamma_x = z0.x*(t-1) + z1.x*t;
        double gamma_y = z0.y*(t-1) + z1.y*t;
        return {gamma_x,gamma_y};
    }

    // Line Integral for Analytic Functions

    template<typename F>
    Complex complex_line_integral_analytic(F f,Complex z0, Complex z1, int N = 2000){
        double delta_t = 1.0 /N; // infinitesimal parameter
        Complex derive_gamma = z1-z0; // derivative of gamma
        Complex sum{0,0};
        for (int k = 0; k<N; ++k) {
            
            double t = k * delta_t; // it actually points of t  
            Complex gamma = complex_curve(z0,z1,t);
            sum = sum + f(gamma) * derive_gamma * delta_t;
        }
        return sum;
    }
}