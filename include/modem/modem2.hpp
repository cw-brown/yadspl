#ifndef MODEM2_HPP
#define MODEM2_HPP

size_t rrcos(size_t k, size_t m, double beta, double* taps){
    size_t n;
    double z, t1, t2, t3, t4, T=1.0;

    double nf, kf, mf;

    double PI = 3.141592;

    size_t h_len = 2 * k * m + 1;

    for(n = 0; n < h_len; ++n){
        nf = (double)n;
        kf = (double)k;
        mf = (double)m;

        z = (nf)/kf - mf;
        t1 = std::cos((1.0+beta)*PI*z);
        t2 = std::sin((1.0-beta)*PI*z);

        if(std::abs(z) < 1e-5){
            taps[n] = 1.0-beta + 4.0*beta/PI; 
        }
        else{
            t3 = 1.0/(4*beta*z);
            double g = 1.0 - 16.0*beta*z*z;
            g*=g;

            if(g < 1e-5){
                double g1, g2, g3, g4;
                g1 = 1.0 + 2.0/PI;
                g2 = std::sin(0.25*PI/beta);
                g3 = 1.0 - 2.0/PI;
                g4 = std::cos(0.25*PI/beta);
                taps[n] = beta/std::sqrt(2.0)*(g1*g2+g3*g4);
            }
            else{
                t4 = 4.0*beta/(PI*std::sqrt(T)*(1.0-(16.0*beta*beta*z*z)));
                taps[n] = t4*(t1 + (t2*t3));
            }
        }

    }
    return h_len;
}

template<typename tap_t, typename out_t>
class fir{
private:
    size_t _n;
    tap_t* _taps;
    size_t _curr;
    out_t* _history;
public:
    constexpr fir(): _n(0), _taps(nullptr), _curr(0), _history(nullptr){}
    constexpr fir(tap_t* taps, const size_t& n): _n(n), _taps(new tap_t[_n]), _curr(0), _history(new out_t[_n]){
        std::uninitialized_copy_n(taps, _n, _taps);
        std::uninitialized_default_construct_n(_history, _n);
    }
    ~fir(){delete[] _taps; delete[] _history;}

    constexpr size_t num_taps() const noexcept{return _n;}
    constexpr tap_t* taps() const noexcept{return _taps;}

    constexpr void set_taps(tap_t* taps, const size_t& n) noexcept{
        if(_taps != nullptr) delete[] _taps;
        if(_history != nullptr) delete[] _history;
        _n = n; _curr = 0;
        _taps = new tap_t[n];
        _history = new out_t[n];
        std::uninitialized_copy_n(taps, _n, _taps);
        std::uninitialized_default_construct_n(_history, _n);
    }
    
    constexpr void insert(const out_t& sample){
        _history[_curr] = sample;
    }
    constexpr out_t operate(){
        out_t accum(0.0);
        size_t idx = _curr;
        for(size_t i = 0; i < _n; ++i){
            accum += _taps[i] * _history[idx];
            idx = idx == 0 ? _n - 1 : idx - 1;
        }
        return accum;
    }
    constexpr void increment(){
        _curr = ++_curr % _n;
    }

    constexpr out_t* filter_n(out_t* input, const size_t& n_in){
        out_t* out = new out_t[n_in];
        for(size_t i = 0; i < n_in; ++i){
            insert(input[i]);
            out[i] = operate();
            increment();
        }
        return out;
    }

    constexpr void reset(){
        std::fill_n(_history, _n, out_t{});
        _curr = 0;
    }
};

template<typename tap_t, typename out_t>
class fir_pfb{
private:
    size_t _prot_size; // total number of filter coefficients
    size_t _arm_size; // total number of coefficients in each arm
    size_t _num_filters;

    out_t* _history;
    size_t _curr;
    tap_t** _taps;
public:
    constexpr fir_pfb():_history(nullptr), _taps(nullptr){}
    constexpr fir_pfb(size_t num_filters, tap_t* taps, size_t n){
        _curr = 0;
        _taps = nullptr;
        _history = nullptr;
        update(num_filters, taps, n);
    }
    ~fir_pfb(){
        delete[] _history;
        for(size_t i = 0; i < _num_filters; ++i) delete[] _taps[i];
        delete[] _taps;
    }

    constexpr size_t prototype_size() const noexcept{return _prot_size;}
    constexpr size_t arm_size() const noexcept{return _arm_size;}
    constexpr size_t num_filters() const noexcept{return _num_filters;}
    constexpr tap_t* arm(size_t idx) const noexcept{return _taps[idx];}

    constexpr void update(size_t num_filters, tap_t* taps, size_t n){
        _num_filters = num_filters;
        _prot_size = n;
        _arm_size = std::ceil(_prot_size / _num_filters);

        if(_taps != nullptr){
            for(size_t i = 0; i < _num_filters; ++i) delete[] _taps[i];
            delete[] _taps;
        }
        _taps = new tap_t*[_num_filters];
        for(size_t i = 0; i < _num_filters; ++i) _taps[i] = new tap_t[_arm_size];
        if(_history != nullptr) delete[] _history;
        _history = new out_t[_arm_size];
        std::uninitialized_default_construct_n(_history, _arm_size);

        tap_t* temp = new tap_t[_num_filters * _arm_size];
        std::uninitialized_default_construct_n(temp, _num_filters * _arm_size);
        std::uninitialized_copy_n(taps, _prot_size, temp);

        for(size_t i = 0; i < _num_filters; ++i){
            for(size_t j = 0; j < _arm_size; ++j){
                _taps[i][j] = temp[i + j * _num_filters];
            }
        }
        reset();
    }

    constexpr void push(out_t sample){
        _history[_curr] = sample;
        _curr = (_curr + 1) % _arm_size;
    }
    constexpr out_t operate(size_t idx){
        out_t accum{};
        size_t index = _curr;
        for(size_t i = 0; i < _arm_size; ++i){
            accum += _taps[idx][i] * _history[index];
            index = index == 0 ? _arm_size - 1 : index - 1;
        }
        return accum;
    }
    constexpr out_t* operate_n(size_t idx, out_t* input, size_t n_in){
        out_t* out = new out_t[n_in];
        for(size_t i = 0; i < n_in; ++i){
            push(input[n_in]);
            out[i] = operate(idx);
        }
        return out;
    }

    constexpr void reset(){
        _curr = 0;
        std::fill_n(_history, _arm_size, out_t{});
    }
};

template<typename base_t, typename out_t>
class fir_resampler{
private:
    size_t P;
    size_t Q;
    base_t _scale;

    fir_pfb<base_t, out_t> _bank;
public:
    fir_resampler(size_t interp, size_t decim, base_t* taps, size_t num_taps){
        P = interp;
        Q = decim;
        _scale = P < Q ? std::sqrt(static_cast<base_t>(P) / static_cast<base_t>(Q)) : 1.0 / std::sqrt(static_cast<base_t>(P) / static_cast<base_t>(Q));
        _bank.update(P, taps, num_taps);
    }

    fir_pfb<base_t, out_t>* bank() noexcept{return &_bank;}

    void operate(out_t* samples, out_t* output){
        // Runs the resampler on a sample array of size Qx1
        // Output placed into array of size Px1
        size_t idx = 0;
        size_t n = 0;
        for(size_t i = 0; i < Q; ++i){
            _bank.push(samples[i]);

            while(idx < P){
                output[n++] = _bank.operate(idx) * _scale;
                idx += Q;
            }
            idx -= P;
        }
    }

    out_t* operate_n(out_t* samples, size_t n_in, size_t& n_written){
        // need to find the closest multiple for Q of the input size
        size_t pad_size = (Q - (n_in % Q)) % Q;
        size_t output_size = pad_size == 0 ? P * n_in / Q : std::floor(n_in * pad_size * P / Q);
        size_t iters = (n_in + pad_size) / Q;
        out_t* out = new out_t[output_size];
        out_t* temp;
        if(pad_size != 0){
            temp = new out_t[n_in + pad_size];
            std::uninitialized_default_construct_n(temp, n_in + pad_size);
            std::copy_n(samples, n_in, temp);
        } 
        else{
            temp = samples;
        }
        for(size_t i = 0; i < iters; ++i){
            operate(temp + i * Q, out + i * P);
        }
        n_written = output_size;
        return out;
    }


};

#endif