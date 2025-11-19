#pragma once
#include <stdexcept>
#include <cmath>
using std::runtime_error;

struct Complex{
    double x; // real part
    double y; // imaginary part

    
    // Addition Operator

    Complex operator+(const Complex& c) const {
        return {x + c.x, y + c.y};
    }

    //Subtraction Operator

    Complex operator-(const Complex &c) const {
        return {x - c.x, y - c.y};
    }

    // Multiplication operator

    Complex operator*(const Complex &c) const {
        return {x*c.x -y*c.y, x*c.y + y*c.x};
    }

    // Multiplication operator for real*complex

    Complex operator*(double s) const {
        return {x * s, y * s};
    }

    // Division operator for complex/real
    
    Complex operator/(double s) const {
    if (s == 0.0) {
        throw std::runtime_error("Division of Complex by zero scalar");
    }
    return {x / s, y / s};
    }

    // Division Operator

    Complex operator/(const Complex &c) const {
        double denom = c.x*c.x + c.y*c.y;

        if(denom ==0.0){
            throw std::runtime_error("Complex division by zero");
        }

        return {(x*c.x+y*c.y)/denom, (y*c.x-x*c.y)/denom};
    }
    // Conjugate of Variable

    Complex operator~() const {
        return{x,-y};
    }

    // Norm Square

    double norm2() const {
        return x*x + y*y;
    }

    // Norm

    double norm() const{
        return std::sqrt(norm2());
    }
};

// Real * Complex (free function)
inline Complex operator*(double s, const Complex& c) {
    return c * s;
}
