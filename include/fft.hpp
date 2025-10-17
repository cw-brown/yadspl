#ifndef RFFT_HPP
#define RFFT_HPP

#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

#include <fftw3.h>

// Reasonably Fast Fourier Transform (in public domain)

// Perform the FFT in place for an array of size 2^k.
// No normalization is done.
void fft_transform_radix2(std::vector<std::complex<double>>& vec, bool inverse);

// Perform the FFT for an arbitrary array using the Bluestein's algorithm.
// The code needs to allocate two supplementary buffers (of size < 4 * n - 3).
// No normalization is done.
void fft_transform_bluestein(std::vector<std::complex<double>>& vec, bool inverse);

// Perform the FFT, choosing the suitable algorithm from the two above.
void fft_transform(std::vector<std::complex<double>>& vec, bool inverse);

#include <complex>

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif 

void fft_transform_radix2(std::vector<std::complex<double>>& vec, bool inverse) {
	size_t n = vec.size();
	int levels = 0;	 // Compute levels = floor(log2(n))
	for (size_t k =	1; (k &	n) == 0; k <<= 1)
		levels++;
	
	// Permute vec by reversing the bits of addresses
	for (size_t i =	0; i < n; i++) {
		// Reverse the bits of i
		size_t j = 0, ii = i;
		for (int k = 0;	k < levels; k++) {
			j = (j << 1) | (ii & 1);
			ii >>= 1;
		}
	
		if (j >	i) {
			std::complex<double> tmp = vec[i];
			vec[i] = vec[j];
			vec[j] = tmp;
		}
	}
	
	// Cooley-Tukey	in place
	for (size_t half = 1; half < n;	half *=	2) {
		size_t size = 2	* half;
		size_t step = n	/ size;
		for (size_t j =	0, k = 0; j < half; j++, k += step) {
			double angle = (inverse	? 2 : -2) * M_PI * k / n;
			std::complex<double> omega = std::polar(1.0, angle);
        		for (size_t i =	0; i < n; i += size) {
				std::complex<double> tmp = vec[i + j + half] * omega;
				vec[i + j + half] = vec[i + j] - tmp;
				vec[i + j] += tmp;
			}
		}
	}
}

void fft_transform_bluestein(std::vector<std::complex<double>>& vec, bool inverse) {
	size_t n = vec.size();
	// Find m = 2^k such that m >= 2 * n + 1
	size_t m = 1;
	while (m <= 2 * n) {
		m *= 2;
	}
	
	// Extended vectors of size m
	std::vector<std::complex<double>> avec(m), bvec(m);
	avec[0]	= vec[0];
	bvec[0]	= 1.0;
	for (size_t i =	1; i < n; i++) {
		size_t k = (i * i) % (2 * n);
		double angle = (inverse	? M_PI : -M_PI)	* k / n;
		std::complex<double> omega = std::polar(1.0, angle);
		avec[i]	= vec[i] * omega;
		bvec[i]	= bvec[m - i] =	std::conj(omega);
	}

	// Convolution
	fft_transform_radix2(avec, false);
	fft_transform_radix2(bvec, false);
	for (size_t i =	0; i < m; i++) {
		avec[i]	*= bvec[i];
	}
	fft_transform_radix2(avec, true);

	for (size_t i =	0; i < n; i++) {
		size_t k = (i * i) % (2 * n);
		double angle = (inverse	? M_PI : -M_PI)	* k / n;
		std::complex<double> omega = std::polar(1.0, angle);
		vec[i] = avec[i] * omega / double(m);
	}
}

void fft_transform(std::vector<std::complex<double>>& vec, bool inverse) {
	size_t n = vec.size();
	if (n <= 1)
		return;
	else if	((n & (n - 1)) == 0)  // Power of 2
		fft_transform_radix2(vec, inverse);
	else
		fft_transform_bluestein(vec, inverse);
}

/**
 * @brief Create a double-sided PSD of a complex FFT transform. This assumes a radix 2 transform
 * 
 * @param vec 
 * @return std::vector<double> 
 */
std::vector<double> make_psd(std::vector<std::complex<double>> vec){
	std::vector<double> output(vec.size());
	size_t half = std::round(vec.size()/2);

	std::vector<double> halfVec(half);

	for(size_t i = 0; i < half; ++i){
		halfVec[i] = 10.0*std::log10(std::abs(vec[i]));
	}
	std::vector<double> halfVecFlip(half);
	std::reverse_copy(halfVec.begin(), halfVec.end(), halfVecFlip.begin());
	auto maximum = *std::max_element(halfVec.begin(), halfVec.end());

	for(size_t i = 0; i < half; ++i){
		output[i] = halfVecFlip[i] - maximum;
	}
	for(size_t i = half; i < vec.size(); ++i){
		output[i] = halfVec[i - half] - maximum;
	}
	return output;
}

double* make_psd(std::complex<double>* vec, size_t n){
	double* output = new double[n];
	size_t half = std::round(n/2);

	std::vector<double> halfVec(half);

	for(size_t i = 0; i < half; ++i){
		halfVec[i] = 10.0*std::log10(std::abs(vec[i]));
	}
	std::vector<double> halfVecFlip(half);
	std::reverse_copy(halfVec.begin(), halfVec.end(), halfVecFlip.begin());
	auto maximum = *std::max_element(halfVec.begin(), halfVec.end());

	for(size_t i = 0; i < half; ++i){
		output[i] = halfVecFlip[i] - maximum;
	}
	for(size_t i = half; i < n; ++i){
		output[i] = halfVec[i - half] - maximum;
	}
	return output;
}

/**
 * @brief Compute the double-sided PSD for a complex signal
 * @param in complex input signal
 * @param n size of the complex input
 * @param output pointer to the magnitude output
 */
double* compute_psd(std::complex<double>* in, size_t n){
	fftw_complex* fft_in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
	fftw_complex* fft_out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));

	for(size_t i = 0; i < n; ++i){
		fft_in[i][0] = in[i].real();
		fft_in[i][1] = in[i].imag();
	}
	fftw_plan plan = fftw_plan_dft_1d(n, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	double* output = new double[n];

	std::complex<double>* fft = reinterpret_cast<std::complex<double>*>(fft_out);
	const size_t half = std::floor(n / 2);
	for(size_t i = 0; i < half; ++i){
		output[i] = 10.0 * std::log10(std::pow(std::norm(fft[half - i]), 2.0) / n);
		// output[i] = 10.0 * std::log10(std::abs(fft[i]) / n);
		output[n - i] = output[i];
	}
	double max_point = *std::max_element(output, output + n);
	for(size_t i = 0; i < n; ++i){
		output[i] -= max_point;
	}

	fftw_destroy_plan(plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
	return output;
}

void real_psd(std::complex<double>* in, size_t n, double* output, double* freqs){
	fftw_complex* fft_in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
	fftw_complex* fft_out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));

	for(size_t i = 0; i < n; ++i){
		fft_in[i][0] = in[i].real();
		fft_in[i][1] = in[i].imag();
	}

	fftw_plan plan = fftw_plan_dft_1d(n, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	std::complex<double>* fft = reinterpret_cast<std::complex<double>*>(fft_out);
	for(size_t i = 0; i < n; ++i){
		output[i] = 10.0 * std::log10(std::abs(fft[i]) / n);
	}
	const size_t half = std::floor(n / 2);
	std::rotate(output, output + half, output + n);
	auto maximum = *std::max_element(output, output + n);

	for(size_t k = 0; k < n; ++k){
		output[k] -= maximum;
        freqs[k] = (k - n / 2.0) / n;
    }
}

#endif // RFFT_IMPLEMENTATION