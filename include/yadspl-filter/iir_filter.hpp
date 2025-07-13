/**
 * @file iir_filter.hpp
 * @author cw-brown (https://github.com/cw-brown)
 * @brief Contains all definitions for building IIR filters, including Biquad sections and bilinear transforms.
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

/**
 * @brief Biquad implements a standard biquad filter with functions to set poles and zeros.
 * @tparam _Tp Floating type of user defined precision.
 */
template<class _Tp>
requires std::is_floating_point_v<_Tp>
class biquad{
public:
    typedef _Tp value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::complex<value_type> complex_type;
    typedef std::pair<complex_type, complex_type> complex_pair;
private:
    static constexpr value_type PI = std::numbers::pi_v<value_type>;
    value_type _a0, _a1, _a2, _b0, _b1, _b2;
    struct state{
        value_type _w1 = value_type{};
        value_type _w2 = value_type{};
        void reset(){_w1 = value_type{}; _w2 = value_type{};}
    } _state;
    struct zpk{
        complex_pair zeros;
        complex_pair poles;
        value_type gain = static_cast<value_type>(1.0);
        zpk(const biquad& s){
            gain = s.getB0()/s.getA0();
            if(s.getA2() == 0 && s.getB2() == 0){
                zeros.first = -s.getB0()/s.getB1();
                zeros.second = 0.0;
                poles.first = -s.getA1();
                poles.second = 0.0;
            } else{
                complex_type c = std::sqrt(complex_type(std::pow(s.getB1(), 2.0) - 4.0*s.getB0()*s.getB2()));
                value_type d = 2.0*s.getB0();
                zeros.first = -(s.getB1() + c)/d;
                zeros.second = (c - s.getB1())/d;

                c = std::sqrt(complex_type(std::pow(s.getA1(), 2.0) - 4.0*s.getA0()*s.getA2()));
                d = 2.0*s.getA0();
                poles.first = -(s.getA1() + c)/d;
                poles.second = (c - s.getA1())/d;
            }
        }
        zpk& operator=(const biquad& s){
            gain = s.getB0()/s.getA0();
            if(s.getA2() == 0 && s.getB2() == 0){
                zeros.first = -s.getB0()/s.getB1();
                zeros.second = 0.0;
                poles.first = -s.getA1();
                poles.second = 0.0;
            } else{
                complex_type c = std::sqrt(complex_type(std::pow(s.getB1(), 2.0) - 4.0*s.getB0()*s.getB2()));
                value_type d = 2.0*s.getB0();
                zeros.first = -(s.getB1() + c)/d;
                zeros.second = (c - s.getB1())/d;

                c = std::sqrt(complex_type(std::pow(s.getA1(), 2.0) - 4.0*s.getA0()*s.getA2()));
                d = 2.0*s.getA0();
                poles.first = -(s.getA1() + c)/d;
                poles.second = (c - s.getA1())/d;
            }
            return *this;
        }
    } _zpk;
public:
    constexpr biquad(): _a0(1.0), _a1(0.0), _a2(0.0), _b0(1.0), _b1(0.0), _b2(0.0), _zpk(*this){}
    constexpr biquad(value_type a0, value_type a1, value_type a2, value_type b0, value_type b1, value_type b2): _a0(a0), _a1(a1/a0), _a2(a2/a0), _b0(b0/a0), _b1(b1/a0), _b2(b2/a0), _zpk(*this){}

    constexpr value_type filter(const value_type& sample){
    // w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
    // y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
        value_type w = sample - _a1*_state._w1 - _a2*_state._w2;
        value_type y = _b0*w + _b1*_state._w1 + _b2*_state._w2;
        _state._w2 = _state._w1;
        _state._w1 = w;
        return y;
    }
    constexpr std::pair<complex_pair, complex_pair> getPoleZeros() const{
        _zpk = *this;
        return std::pair<complex_pair, complex_pair>(_zpk.poles, _zpk.zeros);
    }
    constexpr value_type gain() const{_zpk = *this; return _zpk.gain;}

    constexpr void setCoefficients(value_type a0, value_type a1, value_type a2, value_type b0, value_type b1, value_type b2){
        _a0 = a0;
        _a1 = a1/a0;
        _a2 = a2/a0;
        _b0 = b0/a0;
        _b1 = b1/a0;
        _b2 = b2/a0;
    }

    constexpr void setOnePoleZero(complex_type pole, complex_type zero){
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

    constexpr std::pair<pointer, pointer> response(const size_type num_points){
        std::pair<pointer, pointer> out(new value_type[num_points], new value_type[num_points]);
        value_type d_theta = PI/static_cast<value_type>(num_points);
        for(size_type i = 0; i < num_points; ++i){
            value_type theta = i*d_theta;
            complex_type beta = std::abs(_b0) + std::polar(std::abs(_b1), theta) + std::polar(std::abs(_b2), theta*2.0);
            complex_type alpha = std::abs(_a0) + std::polar(std::abs(_a1), theta) + std::polar(std::abs(_a2), theta*2.0);
            out.first[i] = theta;
            out.second[i] = 20.0*std::log10(std::norm(beta/alpha));
        }
        return out;
    }
};


template<class _Tp>
requires std::is_floating_point_v<_Tp>
class lowpass : public biquad<_Tp>{
public:
    typedef biquad<_Tp>::value_type value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef biquad<_Tp> biquad_type;
public:
    lowpass() = delete;
    lowpass(const value_type& fc, const value_type& q){
        design(fc, q);
    }       

    void design(const value_type& fc, const value_type& q){
        value_type w0 = 2.0*3.141592*fc;
        value_type alpha = std::sin(w0)/2.0/q;
        value_type b0 = (1.0 - std::cos(w0))/2.0;
        value_type b1 = 1.0 - std::cos(w0);
        value_type b2 = b0;
        value_type a0 = 1.0 + alpha;
        value_type a1 = -2.0*std::cos(w0);
        value_type a2 = 1.0 - alpha;
        this->setCoefficients(a0, a1, a2, b0, b1, b2);
    }

};


#endif