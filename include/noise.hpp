#ifndef NOISE_H
#define NOISE_H

#include <algorithm>
#include <random>

/**
 * @brief Generic class for generating different types of signal noise. Intended to work on signals of arithmetic type _T.
 */
template<class _T>
requires std::is_arithmetic_v<_T>
class noise{
private:
    std::mt19937 _gen;
public:
    noise(): _gen(std::random_device{}()){}
    template<class InputIt>
    requires std::input_iterator<InputIt>
    void awgn(InputIt first, InputIt last, const double& snr){
        std::normal_distribution<_T> _d(0.0, 1.0);
        auto wgn = [this, &_d, &snr](_T& x){x += std::sqrt(snr)*_d(_gen);};
        std::for_each(first, last, wgn);
    }
    _T* randomSignal(std::size_t len){
        // Make a signal in the range [-1, 1] of size len
        _T* out = new _T[len];
        std::uniform_real_distribution<_T> dis(-1.0, 1.0);
        for(std::size_t i = 0; i < len; ++i)
            out[i] = dis(_gen);
        return out;
    }
};
#endif