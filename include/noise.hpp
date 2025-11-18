/**
 * @file noise.hpp
 * @author cw-brown (https://github.com/cw-brown)
 * @brief Contains classes dedicated to generating signals and noise types
 * @version 0.1
 * @date 2025-06-18
 */
#ifndef NOISE_H
#define NOISE_H

#include <algorithm>
#include <random>
#include <bitset>
#include "constellations.hpp"

/**
 * @brief Fast noise class for AWGN
 */
class fast_noise{
private:
    std::mt19937 _gen;
    std::normal_distribution<double> _d{0.0, 1.0};
public:
    fast_noise(): _gen(std::random_device{}()){}
    std::complex<double> noise_voltage(double voltage){return voltage * std::complex<double>(_d(_gen), _d(_gen));}

};

class channel_model{
private:
    double _freq_offset;
    double _noise_voltage;
    fast_noise _awgn{};
    std::complex<double> _accum;

    static constexpr double PI = std::numbers::pi;
public:
    channel_model(double frequency_offset, double noise_voltage)
        :_freq_offset(frequency_offset), _noise_voltage(noise_voltage), _accum(1.0){}

    /**
     * @brief Apply offsets and noise to a single sample
     * @param sample 
     */
    void operate(std::complex<double>* sample){
        *sample *= _accum;
        const std::complex<double> w = std::polar(1.0, 2.0 * PI * _freq_offset);
        std::complex<double> noise = _awgn.noise_voltage(_noise_voltage);
        _accum *= w;
        *sample += noise;
    }

    void set_offset(double frequency_offset){_freq_offset = frequency_offset;}
    void set_noise(double noise_power){_noise_voltage = noise_power;}

};

/**
 * @brief Generic class for generating different types of signals. Intended to work on signals of arithmetic type _T.
 */
template<class _T>
requires std::is_arithmetic_v<_T>
class noise{
private:
    std::mt19937 _gen;
public:
    noise(): _gen(std::random_device{}()){}

    /**
     * @brief Creates a random M-ary byte
     * @tparam M 
     * @return std::bitset<M> 
     */
    unsigned int random_m_ary(const int& M){
        std::uniform_int_distribution<unsigned int> dis(0, std::pow(2, M)-1);
        return dis(_gen);
    }

    std::complex<double> random_constellation_point(constellation* constel){
        return constel->get_point(random_m_ary(constel->get_bps()));
    }

    double random_floating(){
        std::normal_distribution<_T> _d(-1.0, 1.0);
        return _d(_gen);
    }

    unsigned int random_bit(){
        std::uniform_int_distribution<unsigned int> dis(0, 1);
        return dis(_gen);
    }

    int random_int_range(int lower, int upper){
        std::uniform_int_distribution<int> dis(lower, upper);
        return dis(_gen);
    }


};
#endif