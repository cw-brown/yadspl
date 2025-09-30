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
     * @brief Corrupts the provided signal with AWGN.
     * @tparam InputIt 
     * @param first 
     * @param last 
     * @param snr 
     */
    template<class InputIt>
    requires std::input_iterator<InputIt>
    void awgn(InputIt first, InputIt last, const double& snr){
        std::normal_distribution<_T> _d(0.0, 1.0);
        auto wgn = [this, &_d, &snr](_T& x){x += std::sqrt(snr)*_d(_gen);};
        std::for_each(first, last, wgn);
    }

    /**
     * @brief Generate random numbers in the range [-1, 1]
     * @param len 
     * @return _T* 
     */
    _T* randomSignal(std::size_t len){
        // Make a signal in the range [-1, 1] of size len
        _T* out = new _T[len];
        std::uniform_real_distribution<_T> dis(-1.0, 1.0);
        for(std::size_t i = 0; i < len; ++i)
            out[i] = dis(_gen);
        return out;
    }

    /**
     * @brief Creates a random M-ary byte
     * @tparam M 
     * @return std::bitset<M> 
     */
    int randomValue(const int& M){
        std::uniform_int_distribution<int> dis(0, std::pow(2, M)-1);
        return dis(_gen);
    }


};
#endif