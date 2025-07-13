/**
 * @file filter.hpp
 * @author cw-brown (https://github.com/cw-brown)
 * @brief Overarching header file for FIR and IIR filters.
 * @version 0.1
 * @date 2025-06-23
 */
#ifndef FILTER_H
#define FILTER_H

#include <algorithm>
#include <memory>
#include <numeric>
#include <numbers>
#include <cmath>
#include <iostream>
#include <complex>
#include <utility>

#include "basic_iterators.hpp"

namespace filter{
typedef double value_type;
typedef std::size_t size_type;
typedef std::ptrdiff_t difference_type;
typedef value_type* pointer;
typedef value_type& reference;
typedef __iterator<value_type> iterator;
typedef __const_iterator<value_type> const_iterator;
typedef __reverse_iterator<value_type> reverse_iterator;
typedef __const_reverse_iterator<value_type> const_reverse_iterator;
typedef std::complex<value_type> complex_type;

/**
 * @brief Struct for containing a filter response.
 */
struct filter_response{
    pointer frequencies;
    pointer magnitude;
    pointer phase;
};

enum class filter_type{
    Lowpass,
    Highpass,
    Bandpass,
    Bandstop
};

enum class window{
    None,
    Triangular,
    Von_Hann,
    Hamming,
    Blackman_Harris
};

static constexpr value_type PI = std::numbers::pi_v<value_type>;
static constexpr value_type TWOPI = 2.0*PI;
static constexpr value_type EPS = 1e-6;

template<filter_type _F, window _W = window::None>
class fir{
private:
    size_type _n;
    pointer _taps;

    void synthesize() requires(_W == window::None){}

    void synthesize() requires(_W == window::Triangular){
        value_type mid = static_cast<value_type>(_n - 1)/2.0;
        for(size_type i = 0; i < _n; ++i) _taps[i] *= 1.0 - std::abs((i - mid)/mid);
    }

    void synthesize() requires(_W == window::Von_Hann){
        for(size_type i = 0; i < _n; ++i) _taps[i] *= 0.5*(1.0 - std::cos((2.0*PI*i)/(_n - 1)));
    }

    void synthesize() requires(_W == window::Hamming){
        for(size_type i = 0; i < _n; ++i) _taps[i] *= 0.54 - 0.46*std::cos((2.0*PI*i)/(_n - 1));
    }

    void synthesize() requires(_W == window::Blackman_Harris){
        for(size_type n = 0; n < _n; ++n) _taps[n] *= 0.35875 - 0.48829*std::cos(2.0*PI*n/_n) + 0.14128*std::cos(4.0*PI*n/_n) - 0.01168*std::cos(6.0*PI*n/_n);
    }

public:
    fir(): _n(0), _taps(nullptr){}
    fir(fir& other): _n(other._n), _taps(new value_type[other._n]){std::uninitialized_copy(other.begin(), other.end(), begin());}
    fir(fir&& other): _n(other._n), _taps(std::move(other._taps)){other._taps = nullptr;}

    void design(size_type order, value_type fc)
    requires(_F == filter_type::Lowpass){
        _n = order;
        _taps = new value_type[_n];
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = fc / PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = std::sin(m*fc)/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
        synthesize();
    }

    void design(size_type order, value_type fc)
    requires(_F == filter_type::Highpass){
        _n = order;
        _taps = new value_type[_n];
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = 1.0 - fc / PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = -std::sin(m*fc)/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
        synthesize();
    }

    void design(size_type order, value_type fc1, value_type fc2)
    requires(_F == filter_type::Bandpass){
        _n = order;
        _taps = new value_type[_n];
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = (fc2 - fc1)/PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = (std::sin(m*fc2) - std::sin(m*fc1))/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
        synthesize();
    }

    void design(size_type order, value_type fc1, value_type fc2)
    requires(_F == filter_type::Bandstop){
        _n = order;
        _taps = new value_type[_n];
        difference_type n_max = order % 2 ? (order - 1)/2 : order/2;
        if(order % 2) _taps[n_max] = 1.0 + (fc1 - fc2)/PI;
        for(difference_type n = 0; n < n_max; ++n){
            difference_type m = n - n_max;
            _taps[n] = (std::sin(m*fc1) - std::sin(m*fc2))/m/PI;
            _taps[order - 1 - n] = _taps[n];
        }
        synthesize();
    }

    filter_response response(size_type freqs){
        filter_response out;
        out.frequencies = new value_type[freqs];
        out.magnitude = new value_type[freqs];
        out.phase = new value_type[freqs];
        value_type d_theta = PI/static_cast<value_type>(freqs);
        for(size_type i = 0; i < freqs; ++i){
            value_type theta = i*d_theta;
            std::complex<value_type> work(0.0, 0.0);
            for(size_type j = 0; j < _n; ++j){
                work += std::polar(_taps[j], 1.0*j*theta);
            }
            out.frequencies[i] = theta;
            out.magnitude[i] = 20.0*std::log10(std::norm(work));
            out.phase[i] = std::arg(work);
        }
        return out;
    }

    pointer filter_n(pointer input, size_type m){
        pointer out = new value_type[m];
        for(size_type i = 0; i < m; ++i){
            value_type work = value_type{};
            for(size_type j = 0; j < _n; ++j){
                if(i >= j) work += _taps[j] * input[i - j];
            }
            out[i] = work;
        }
        return out;
    }

    iterator begin(){return iterator(_taps[0]);}
    const_iterator begin() const{return const_iterator(_taps[0]);}
    iterator end(){return iterator(_taps[_n]);}
    const_iterator end() const{return const_iterator(_taps[_n]);}
};

namespace iir{

class biquad{
private:
    value_type _a0, _a1, _a2, _b0, _b1, _b2;
    struct state{
        value_type _w1{};
        value_type _w2{};
        void reset(){_w1 = value_type{}; _w2 = value_type{};}
    } _state;
public:
    biquad(): _a0(1.0), _a1(0.0), _a2(0.0), _b0(1.0), _b1(0.0), _b2(0.0){}
    biquad(value_type a0, value_type a1, value_type a2, value_type b0, value_type b1, value_type b2): _a0(a0), _a1(a1/a0), _a2(a2/a0), _b0(b0/a0), _b1(b1/a0), _b2(b2/a0){}

    value_type filter(value_type& sample){
        value_type w = sample - _a1*_state._w1 - _a2*_state._w2;
        value_type y = _b0*w + _b1*_state._w1 + _b2*_state._w2;
        _state._w2 = _state._w1;
        _state._w1 = w;
        return y;
    }

    void setCoefficients(value_type a0, value_type a1, value_type a2, value_type b0, value_type b1, value_type b2){
        _a0 = a0;
        _a1 = a1/a0;
        _a2 = a2/a0;
        _b0 = b0/a0;
        _b1 = b1/a0;
        _b2 = b2/a0;
    }

    void setOnePoleZero(complex_type pole, complex_type zero){
        if(pole.imag() != 0.0) throw std::invalid_argument("biquad::setOnePoleZero: pole imaginary component is non-zero");
        if(zero.imag() != 0.0) throw std::invalid_argument("biquad::setOnePoleZero: zero imaginary component is non-zero");
        setCoefficients(1.0, -pole.real(), 0.0, 1.0, -zero.real(), 0.0);
    }

    constexpr void setAllPoleZeros(complex_type p1, complex_type p2, complex_type z1, complex_type z2){
        value_type a1, a2, b1, b2;
        if(p1.imag() != 0){
            if(p2 != std::conj(p1)) throw std::invalid_argument("biquad::setAllPoleZeros: poles must be complex conjugate pairs");
            a1 = -2.0*p1.real();
            a2 = std::norm(p1);
        } else{
            if(p2.imag() != 0) throw std::invalid_argument("biquad::setAllPoleZeros: poles must be complex conjugate pairs");
            a1 = -(p1.real() + p2.real());
            a2 = p1.real()*p2.real();
        }
        if(z1.imag() != 0){
            if(z2 != std::conj(z1)) throw std::invalid_argument("biquad::setAllPoleZeros: zeros must be complex conjugate pairs");
            b1 = -2.0*z1.real();
            b2 = std::norm(z1);
        } else{
            if(z2.imag() != 0) throw std::invalid_argument("biquad::setAllPoleZeros: zeros must be complex conjugate pairs");
            b1 = -(z1.real() + z2.real());
            b2 = z1.real()*z2.real();
        }
        setCoefficients(1.0, a1, a2, 1.0, b1, b2);
    }

    constexpr value_type getA0() const{return _a0;}
    constexpr value_type getA1() const{return _a1*_a0;}
    constexpr value_type getA2() const{return _a2*_a0;}
    constexpr value_type getB0() const{return _b0*_a0;}
    constexpr value_type getB1() const{return _b1*_a0;}
    constexpr value_type getB2() const{return _b2*_a0;}

    filter_response response(size_type freqs){
        filter_response out;
        out.frequencies = new value_type[freqs];
        out.magnitude = new value_type[freqs];
        out.phase = new value_type[freqs];
        value_type d_theta = PI/static_cast<value_type>(freqs);
        for(size_type i = 0; i < freqs; ++i){
            value_type theta = i*d_theta;
            complex_type beta = _b0 + std::polar(std::abs(_b1), theta) + std::polar(std::abs(_b2), 2.0*theta);
            complex_type alpha = _a0 + std::polar(std::abs(_a1), theta) + std::polar(std::abs(_a2), 2.0*theta);
            out.frequencies[i] = theta;
            out.magnitude[i] = 20.0*std::log10(std::norm(beta/alpha));
            out.phase[i] = std::arg(beta/alpha);
        }
        return out;
    }
};

template<size_type order, filter_type _F>
class butterworth{
private:
    biquad _2nd_order_biquad;
    biquad _1st_order_biquad;
public:
    butterworth(value_type fc, value_type fs){
        /**
         * Odd order filter: design a 1st and 2nd order butterworth prototype, use bilinear transform, the filter will recursively use these to filter samples
         * Even order filter: design a 2nd order butterworth prototype
         */

        if(order % 2){
            // odd filter
            value_type T = 1.0/fs;
            
        }
    }
};

}
}
#endif