#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <typeinfo>
#include <memory>
#include <cassert>

// #include "implot.h"

// #include "helper_funcs.h"
// #include "fir_filter.hpp"
#include "polynomial.hpp"

int main(){
    poly<double, 3> a = {-12.0, 15.0, 1.0, -2.0}; // -12 + 15x + 1x^2 - 2x^3
    poly<double, 1> b = {7.0, 2.0}; // 7 + 2x
    poly<double, 0> c = {2.0}; // 2.0

    assert
    (""
        "Addition of polynomials:" &&
        (a + b) == (poly<double, 3>{-5.0, 17.0, 1.0, -2.0}) && 
        (b + a) == (poly<double, 3>{-5.0, 17.0, 1.0, -2.0}) && 
        (b + b) == (poly<double, 1>{14.0, 4.0}) &&
        (c + c) == (poly<double, 0>{4.0}) &&
    "");




    return 0;
}