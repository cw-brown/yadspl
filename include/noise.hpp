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
    fast_noise _awgn{};


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