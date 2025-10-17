#include <iostream>
#include <complex>
#include <array>
#include <memory>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <utility>

#include <fftw3.h>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
#include "noise.hpp"
#include "yadpsl_math.hpp"
#include "constellations.hpp"
#include "symbol_rec.hpp"
#include "ring.hpp"
#include "polyphase.hpp"
#include "graphics.hpp"
#include "modem.hpp"

#include "fft.hpp"

#define DO_WINDOW true

int main(){
    size_t sps = 4;
    size_t n_filt = 32;
    size_t points = 250;
    size_t N = sps * points;
    size_t k = 0;

    static float v = 0.0f;
    static float offset = 0.0f;

    constellation_qpsk QPSK{};
    noise<double> sig_gen{};
    static channel_model model(offset, v);
    rectangular_modulator modulator(&QPSK, sps, n_filt, 0.35);
    std::complex<double>* sig_buffer = new std::complex<double>[N];
    double* sig_psd = new double[N];
    double* sig_freqs = new double[N];
    std::complex<double>* sig_buff_copy = new std::complex<double>[N];

    firefighter recovery(&QPSK, sps, n_filt, 2.0*3.14/200.0, 0.35);
    std::complex<double>* rec_buffer = new std::complex<double>[N];
    double* rec_psd = new double[N];
    double* rec_freqs = new double[N];

    for(size_t i = 0; i < points; ++i){
        auto point = sig_gen.random_int_range(0, QPSK.get_size() - 1);
        modulator.operate(point, sig_buffer + (i * sps));
    }

    auto prot = root_nyquist(n_filt, n_filt * sps, 1.0, 0.35, 8 * n_filt * sps);
    resampler<std::complex<double>> decim(1.0 / sps, n_filt, prot);
    std::complex<double>* symb_buffer = new std::complex<double>[points];
    std::complex<double>* symb = new std::complex<double>[sps];
    double* symb_real = new double[points];
    double* symb_imag = new double[points];

#if DO_WINDOW
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        std::copy_n(sig_buffer, N, sig_buff_copy);
        for(size_t i = 0; i < N; ++i){
            model.operate(sig_buff_copy + i);
        }
        real_psd(sig_buff_copy, N, sig_psd, sig_freqs);

        for(size_t i = 0; i < N; ++i){
            recovery.operate(sig_buff_copy[i], rec_buffer + i);
        }
        real_psd(rec_buffer, N, rec_psd, rec_freqs);

        int j = 0;
        for(size_t i = 0; i < N; ++i){
            int n = decim.operate(rec_buffer[i], symb);
            if(n > 0){
                symb_buffer[j] = *symb;
                symb_real[j] = symb->real();
                symb_imag[j] = symb->imag();
                j++;
            }
        }

   

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        if(ImGui::BeginTabBar("Main Tabs")){
        if(ImGui::BeginTabItem("Item 0")){
            ImGui::SliderFloat("Offset", &offset, 0.0, 1.0);
            ImGui::SliderFloat("Noise", &v, 0.0, 1.0);
            model.set_noise(v);
            model.set_offset(offset);
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                ImPlot::SetupAxesLimits(-0.5, 0.5, -50, 5);
                ImPlot::PlotLine("Signal FFT", sig_freqs, sig_psd, N);
                ImPlot::PlotLine("FLL FFT", rec_freqs, rec_psd, N);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Plot 1", ImVec2(-1, 750))){
                ImPlot::PlotScatter("", symb_real, symb_imag, points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 1")){
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 2")){
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        }
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();
#endif
    return 0;
}

