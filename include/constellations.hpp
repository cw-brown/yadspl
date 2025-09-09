#ifndef CONSTELLATION_H
#define CONSTELLATION_H

#include <vector>
#include <complex>

class constellation{
protected:
    std::vector<std::complex<double>> _points;
    unsigned int _n;

    double distance(unsigned int idx, const std::complex<double>& sample) const{
        return std::norm(sample - _points.at(idx));
    }
public:
    enum normalization{
        power,
        amplitude
    };

    constellation(const std::vector<std::complex<double>>& points)
        : _points(points), _n(points.size()){}
    constellation()
        : _points(1), _n(0){}

    virtual ~constellation(){}

    std::complex<double> get_point(const unsigned int& value) const{return _points.at(value);}
    unsigned int get_size() const{return _n;}
    std::vector<std::complex<double>> get_constellation() const{return _points;}
    unsigned int get_bps() const{return std::log2(_n);}
    std::pair<double*, double*> get_IQ(){
        double* I = new double[_n];
        double* Q = new double[_n];
        for(size_t i = 0; i < _n; ++i){
            I[i] = _points.at(i).real();
            Q[i] = _points.at(i).imag();
        }
        return std::pair<double*, double*>(I, Q);
    }

    // virtual unsigned int decision(const std::complex<double>& sample) = 0;

    void normalize(normalization norm){
        switch(norm){
        case power:{
            double pow = 0;
            for(auto&& v : _points) pow += std::norm(v);
            double scale = std::sqrt(pow / _n);
            for(auto&& v : _points) v /= scale;
            break;
        }
        case amplitude:{
            double amp = 0;
            for(auto&& v : _points) amp += std::abs(v);
            double scale = _n / amp;
            for(auto&& v : _points) v *= scale;
            break;
        }
        }
    }
    unsigned int closest_point(const std::complex<double>& sample){
        double cur_min = distance(0, sample);
        for(unsigned int i = 1; i < _n; ++i){
            double dist = distance(i, sample);
            cur_min = dist < cur_min ? dist : cur_min;
        }
        return cur_min;
    }

};



/*
1000   1101 | 1100   1001
           |
1111   1010 | 1011   1110
   -----------------
0100   0001 | 0000   0101
           |
0011   0110 | 0111   0010
*/
class constellation_16qam : public constellation{
public:
    constellation_16qam(){
        _points.resize(16);
        _points[0] = std::complex<double>(1, -1);
        _points[1] = std::complex<double>(-1, -1);
        _points[2] = std::complex<double>(3, -3);
        _points[3] = std::complex<double>(-3, -3);
        _points[4] = std::complex<double>(-3, -1);
        _points[5] = std::complex<double>(3, -1);
        _points[6] = std::complex<double>(-1, -3);
        _points[7] = std::complex<double>(1, -3);
        _points[8] = std::complex<double>(-3, 3);
        _points[9] = std::complex<double>(3, 3);
        _points[10]= std::complex<double>(-1, 1);
        _points[11]= std::complex<double>(1, 1);
        _points[12]= std::complex<double>(1, 3);
        _points[13]= std::complex<double>(-1, 3);
        _points[14]= std::complex<double>(3, 1);
        _points[15]= std::complex<double>(-3, 1);
        _n = 16;
        normalize(normalization::power);
    }

    ~constellation_16qam() override{}

    // unsigned int decision(const std::complex<double>& sample) override{
    //     return 1;
    // }

};

/*
01 | 11
-------
00 | 10
*/
class constellation_qpsk : public constellation{
public:
    constellation_qpsk(){
        _points.resize(4);
        _points[0] = std::complex<double>(-1, -1);
        _points[1] = std::complex<double>(-1, 1);
        _points[2] = std::complex<double>(1, -1);
        _points[3] = std::complex<double>(1, 1);
        _n = 4;
        normalize(normalization::power);
    }

    ~constellation_qpsk(){}
};

class constellation_bpsk : public constellation{
public:
    constellation_bpsk(){
        _points.resize(2);
        _points[0] = std::complex<double>(1, 0);
        _points[1] = std::complex<double>(-1, 0);
        _n = 2;
        normalize(normalization::power);
    }

    ~constellation_bpsk(){}
};



#endif