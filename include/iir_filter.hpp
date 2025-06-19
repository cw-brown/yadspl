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
    Butterworth_Lowpass,
    None
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
private:
    size_type _nz;
    size_type _mp;
    std::complex<value_type>* _zeros;
    std::complex<value_type>* _poles;
    static constexpr value_type PI = std::numbers::pi_v<__T>;
    static constexpr value_type eps = static_cast<__T>(1e-5);    
public:
    constexpr iir_filter(): _nz(0), _mp(0), _zeros(nullptr), _poles(nullptr){}
    constexpr iir_filter(const size_type& order, const value_type& fc) requires(__F == iir_type::Butterworth_Lowpass){
        // Design an iir filter from a prototype butterworth using the bilinear transformation
        // all fc are mormalized 0-pi

        using namespace std::literals::complex_literals;
        // Prototype
        _nz = 0;
        _mp = order;
        _zeros = new std::complex<value_type>[1];
        _poles = new std::complex<value_type>[order];

        _zeros[0] = std::pow(fc, 2.0);
        for(size_type i = 0; i < order; ++i){
            _poles[i] = fc*std::exp(1.0i*PI*(2.0*(i + 1) + static_cast<value_type>(order) - 1.0)/2.0/static_cast<value_type>(order));
        }
    }

    constexpr std::complex<value_type>* getPoles() const noexcept{return _poles;}
    constexpr std::complex<value_type>* getZeros() const noexcept{return _zeros;}
    constexpr size_type getNumPoles() const noexcept{return _mp;}
    constexpr size_type getNumZeros() const noexcept{return _nz;}

    constexpr std::pair<pointer, pointer> response(const size_type& points){
        using namespace std::literals::complex_literals;
        std::pair<pointer, pointer> out(new value_type[points], new value_type[points]);
        value_type d_theta = PI/static_cast<value_type>(points);
        for(size_type idx = 0; idx < points; ++idx){
            value_type theta = idx*d_theta;
            std::complex<value_type> numer(1.0, 0.0);
            std::complex<value_type> denom(1.0, 0.0);
            if(_nz == 0) numer = _zeros[0] + 0.0i;
            for(size_type n = 0; n < _nz; ++n){
                numer *= std::exp(static_cast<value_type>(n)*1.0i*theta) + _zeros[n];
            }
            for(size_type m = 0; m < _mp; ++m){
                denom *= std::exp(static_cast<value_type>(m)*1.0i*theta) + _poles[m];
            }
            out.first[idx] = theta;
            out.second[idx] = 20.0 * std::log10(std::norm(numer/denom));
        }
        return out;
    }

};
#endif