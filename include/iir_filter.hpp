/**
 * @file iir_filter.hpp
 * @author cw-brown (https://github.com/cw-brown)
 * @brief IIR filter class for some basic filter types.
 * @version 0.1
 * @date 2025-06-18
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef IIR_FILT_H
#define IIR_FILT_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>
#include <numeric>
#include <memory>

#include "basic_iterators.hpp"

enum class iir_type{
    Chebyshev,
    Highpass
};

template<class __T, iir_type __F>
requires std::is_floating_point_v<__T>
class iir_filter{
public:
    typedef __T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef __T* pointer;
    typedef __T& reference;
    typedef __iterator<__T> iterator;
    typedef __const_iterator<__T> const_iterator;
    typedef __reverse_iterator<__T> reverse_iterator;
    typedef __const_reverse_iterator<__T> const_reverse_iterator;
    typedef std::allocator_traits<std::allocator<__T>>::allocator_type allocator_type;
private:
    size_type _n;
    size_type _m;
    pointer _a;
    pointer _b;
    std::allocator<__T> _alloc;
    static constexpr value_type PI = std::numbers::pi_v<__T>;
    static constexpr value_type eps = static_cast<__T>(1e-5);
public:
    
};
#endif