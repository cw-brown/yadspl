#include <cassert>

#include "polynomial.hpp"

int main(){
    poly<double, 3> a = {-12.0, 15.0, 1.0, -2.0}; // -12 + 15x + 1x^2 - 2x^3
    poly<double, 1> b = {7.0, 2.0}; // 7 + 2x
    poly<double, 0> c = {2.0}; // 2.0

    assert
    (""
        "Addition operations:" &&
        (a + b) == (poly<double, 3>{-5.0, 17.0, 1.0, -2.0}) && 
        (b + a) == (poly<double, 3>{-5.0, 17.0, 1.0, -2.0}) && 
        (b + b) == (poly<double, 1>{14.0, 4.0}) &&
        (b + c) == (poly<double, 1>{9.0, 2.0}) &&
        (c + c) == (poly<double, 0>{4.0}) &&
        (b + 2.0) == (poly<double, 1>{9.0, 2.0}) &&
        (c + 2.0) == (poly<double, 0>{4.0}) &&
        (2.0 + b) == (poly<double, 1>{9.0, 2.0}) &&
        (2.0 + c) == (poly<double, 0>{4.0}) &&
    "");

    assert
    (""
        "Subtraction operations:" &&
        (a - b) == (poly<double, 3>{-19.0, 13.0, 1.0, -2.0}) &&
        (b - a) == (poly<double, 3>{19.0, -13.0, -1.0, 2.0}) &&
        (b - b) == (poly<double, 1>{0.0, 0.0}) &&
        (b - c) == (poly<double, 1>{5.0, 2.0}) &&
        (c - c) == (poly<double, 0>{0.0}) && 
        (b - 2.0) == (poly<double, 1>{5.0, 2.0}) &&
        (c - 2.0) == (poly<double, 0>{0.0}) &&
        (2.0 - b) == (poly<double, 1>{-5.0, -2.0}) &&
        (2.0 - c) == ((poly<double, 0>{0.0})) &&
    "");

    assert
    (""
        "Multiplication operations:" &&
        (a * b) == (poly<double, 4>{-84.0, 81.0, 37.0, -12.0, -4.0}) &&
        (b * a) == (poly<double, 4>{-84.0, 81.0, 37.0, -12.0, -4.0}) &&
        (b * b) == (poly<double, 2>{49.0, 28.0, 4.0}) &&
        (b * c) == (poly<double, 1>{14.0, 4.0}) &&
        (c * c) == (poly<double, 0>{4.0}) &&
        (b * 2.0) == (poly<double, 1>{14.0, 4.0}) &&
        (c * 2.0) == (poly<double, 0>{4.0}) &&
        (2.0 * b) == (poly<double, 1>{14.0, 4.0}) &&
        (2.0 * c) == (poly<double, 0>{4.0}) &&
    "");

    assert
    (""
        "Division Operations:" &&
        (b / b) == (poly<double, 1>{1.0, 1.0}) &&
        (b / 2.0) == (poly<double, 1>{7.0/2.0, 1.0}) &&
    "");
}