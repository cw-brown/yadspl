#ifndef POLYPHASE_H
#define POLYPHASE_H

#include <vector>
#include <complex>
#include <numbers>
#include <numeric>
#include <cmath>

/**
 * @brief Returns the taps for a root nyquist filter.
 * @param gain total gain
 * @param fs sampling frequency
 * @param sr symbol rate (symbols/s)
 * @param alpha excess bandwidth/rolloff 
 * @param n number of taps
 * @return std::vector<double> 
 */
std::vector<double> root_nyquist(const double& gain, const double& fs, const double& sr, const double& alpha, const size_t& n){
    if(gain <= 0.0) throw std::domain_error("make_nyquist: gain cannot be less than 0");
    if(alpha <= 0.0 || alpha > 1.0) throw std::domain_error("make_nyquist: rolloff must be in range (0, 1]");

    const size_t n_taps = n | 1;
    std::vector<double> taps(n_taps, 0.0);
    double power = 0.0;
    const double sps = fs/sr;
    const double PI = std::numbers::pi;

    for(size_t i = 0; i < n_taps; ++i){
        double num, denom;
        double xidx = (double)i - (double)n_taps/2.0;
        double x1 = PI*xidx/(double)sps;
        double x2 = 4.0*alpha*xidx/(double)sps;
        double x3 = std::pow(x2, 2.0) - 1.0;

        if(std::abs(x3) >= 1.0e-5){
            if(i != n_taps/2){
                num = std::cos((1.0 + alpha)*x1) 
                    + std::sin((1.0 - alpha)*x1)/(4.0*alpha*xidx/(double)sps);
            }
            else{
                num = std::cos((1.0 + alpha)*x1)
                    + (1.0 - alpha)*PI/(4.0*alpha);
            }
            denom = x3*PI;
        }
        else{
            if(alpha == 1.0){
                taps[i] = -1.0;
                power += taps[i];
                continue;
            }
            x3 = (1.0 - alpha)*x1;
            x2 = (1.0 + alpha)*x1;
            num = std::sin(x2)*(1.0 + alpha)*PI
                - std::cos(x3)*((1.0 - alpha)*PI*(double)sps)/(4.0*alpha*xidx)
                + std::sin(x3)*std::pow((double)sps, 2.0)/(4.0*alpha*std::pow(xidx, 2.0));
            denom = -32.0*PI*std::pow(alpha, 2.0)*xidx/(double)sps;
        }
        taps[i] = 4.0*alpha*num/denom;
        power += taps[i];
    }
    for(size_t i = 0; i < n_taps; ++i){
        taps[i] = taps[i]*gain/power;
    }
    return taps;
}

/**
 * @brief Container class for taps of an FIR filter.
 * 
 */
class fir_taps{
private:
    std::vector<double> _taps;
    size_t _n;
    std::vector<std::complex<double>> _history;
    size_t _curr = 0;
public:
    fir_taps(const std::vector<double>& taps, const size_t& n): _taps(taps), _n(n), _history(n, 0.0){}
    std::complex<double> filter1(const std::complex<double>& sample){
        std::complex<double> accum(0.0, 0.0);
        _history[_curr] = sample;
        size_t idx = _curr;
        for(size_t i = 0; i < _n; ++i){
            accum += _taps[i]*_history[idx];
            if(idx == 0) idx = _history.size() - 1;
            else idx -= 1;
        }
        _curr = (_curr + 1)%_history.size();
        return accum;
    }
    void update_taps(const std::vector<double>& taps){_taps = taps; _n = taps.size();}

    void reset(){
        _history.assign(_n, 0.0);
        _curr = 0;
    }

    std::vector<double> taps() const{return _taps;}
    size_t size() const{return _n;}
};

/**
 * @brief polyphase_filter_bank is an internal class of symbol recovery that peforms as a polyphase filter bank. It is made with a prototype, 
 * usually a root nyquist filter.
 */
class polyphase_filter_bank{
private:
    size_t _n; // number of arms of the bank
    size_t _taps_per_filter;
    std::vector<fir_taps> _filters;
    std::vector<fir_taps> _deriv_filters;
    std::vector<double> _prototype;
    std::vector<double> _deriv_prototype;

    void set_taps(){
        std::vector<double> prot_filled(_prototype);
        while(prot_filled.size() < _n * _taps_per_filter) prot_filled.push_back(0.0);
        for(size_t i = 0; i < _n; ++i){
            std::vector<double> tmp(_taps_per_filter, 0.0);
            for(size_t j = 0; j < _taps_per_filter; ++j){
                tmp[j] = prot_filled[i + j*_n];
            }
            _filters[i].update_taps(tmp);
        }

        std::vector<double> deriv_filled(_deriv_prototype);
        while(deriv_filled.size() < _n * _taps_per_filter) deriv_filled.push_back(0.0);
        for(size_t i = 0; i < _n; ++i){
            std::vector<double> tmp(_taps_per_filter, 0.0);
            for(size_t j = 0; j < _taps_per_filter; ++j){
                tmp[j] = deriv_filled[i + j*_n];
            }
            _deriv_filters[i].update_taps(tmp);
        }
    };

    void make_derivative(){
        _deriv_prototype.reserve(_prototype.size());
        _deriv_prototype.push_back(0.0);
        std::vector<double> diff{-1.0, 0.0, 1.0};
        double power = 0.0;
        for(size_t i = 0; i < _prototype.size() - 2; ++i){
            double accum = 0.0;
            for(size_t j = 0; j < diff.size(); ++j){
                accum += diff[j]*_prototype[i + j];
            }
            _deriv_prototype.push_back(accum);
            power += std::abs(accum);
        }
        _deriv_prototype.push_back(0.0);
        for(size_t i = 0; i < _deriv_prototype.size(); ++i){
            _deriv_prototype[i] *= _n/power;
        }
    };
public:
    polyphase_filter_bank(): _n(0), _taps_per_filter(0){}
    constexpr polyphase_filter_bank(const std::vector<double>& prototype, const size_t& arms){
        update(prototype, arms);
    }

    size_t get_size() const{return _n;}
    size_t get_taps_per_arm() const{return _taps_per_filter;}
    fir_taps get_arm(const size_t& idx) const{return _filters.at(idx);}
    fir_taps get_deriv_arm(const size_t& idx) const{return _deriv_filters.at(idx);}
    std::vector<double> get_prototype() const{return _prototype;}
    std::vector<double> get_deriv_prototype() const{return _deriv_prototype;}

    std::complex<double> filter(const std::complex<double>& sample, const size_t& arm){
        return _filters[arm].filter1(sample);
    }
    std::complex<double> deriv_filter(const std::complex<double>& sample, const size_t& arm){
        return _deriv_filters[arm].filter1(sample);
    }

    constexpr void update(const std::vector<double>& prototype, const size_t& arms){
        _n = arms;
        _prototype = prototype;
        make_derivative();

        // Set up our tap vectors
        _taps_per_filter = std::ceil(prototype.size()/_n);
        _filters.reserve(_n);
        _deriv_filters.reserve(_n);
        fir_taps tmp_taps(std::vector<double>(1, 0.0), 1);
        for(size_t i = 0; i < _n; ++i){
            _filters.emplace_back(tmp_taps);
            _deriv_filters.emplace_back(tmp_taps);
        }

        set_taps();
    }

    void clear(){
        _filters.clear();
        _deriv_filters.clear();
        _prototype.clear();
        _deriv_prototype.clear();
    }
};

/**
 * @brief Polyphase upsampler is an FIR polyphase bank up sampling filter for digital modulation.
 * 
 */
class polyphase_upsampler{
private:
    size_t _up_rate; // interpolation rate
    size_t _down_rate; // decimation rate
    double _fractional_rate; // fractional rate
    size_t _curr; // current arm
    size_t _n_filts; // number of filters in the bank
    double _accum;

    polyphase_filter_bank _bank;
public:
    polyphase_upsampler(const double& rate, const std::vector<double>& prototype, size_t nfilts)
        : _n_filts(nfilts), _bank(prototype, nfilts){
        // Set the different rates for the filter
        _up_rate = nfilts;
        _down_rate = std::floor(_up_rate / rate);
        _fractional_rate = _up_rate / rate - _down_rate;

        _curr = std::floor(_n_filts / 2);
        _accum = 0.0;
    }

    /**
     * @brief Filter N samples 
     * 
     * @tparam inputIt 
     * @param it Iterator start for the block of samples
     * @param n Number of samples to filter
     */
    std::vector<std::complex<double>> filterN(std::complex<double>* samples, size_t n){
        std::vector<std::complex<double>> output(n * _up_rate);
        
        // for(size_t i = 0; i < n; ++i){
        //     while(_curr < _up_rate){
        //         std::complex<double> x1 = _bank.filter(samples[i])
        //     }
        
        // }

        size_t i = 0; 
        size_t j = 0;
        while(i < n){
            while(_curr < _up_rate){
                std::complex<double> x1 = _bank.filter(samples[i], _curr);
                std::complex<double> x2 = _bank.deriv_filter(samples[i], _curr);

                output[j] = x1 + x2 * _accum;
                j++;
                _accum += _fractional_rate;
                _curr += _down_rate + std::floor(_accum);
                _accum = std::fmod(_accum, 1.0);
            }
            i += _curr / _up_rate;
            _curr %= _up_rate;
        }
        return output;
    }
};

#endif