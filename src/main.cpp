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
    using cpx = std::complex<double>;
    size_t sps = 2;
    size_t n_filt = 32;
    size_t points = 250;
    size_t preamble_len = 0;
    size_t N = preamble_len + sps * points;

    static float v = 0.0f;
    static float offset = 0.0f;

    // static float alpha = 0.71419;
    // static float beta = 0.83998;
    // static float max_freq = 0.13;

    constellation_bpsk constel{};
    noise<double> sig_gen{};
    channel_model channel(offset, v);
    rectangular_modulator modulator(&constel, sps, n_filt, 0.35);
    firefighter recovery(&constel, sps, n_filt, 2.0*3.1415/100.0, 0.35);
    cpx* data = new cpx[N];
    cpx* data_c = new cpx[N];

    unsigned int* pattern = new unsigned int[points];

    static int symb_delay = 12;

    for(size_t i = 0; i < preamble_len; ++i){
        int val;
        if(i % 2){
            val = 0;
        }
        else{
            val = 1;
        }
        modulator.operate(val, data + i * sps);
        pattern[i] = val;

    }

    for(unsigned int i = preamble_len; i < points; ++i){
        // auto point = (i * i * 2) % (constel.get_size());
        auto point = sig_gen.random_int_range(0, constel.get_size() - 1);
        modulator.operate(point, data + i * sps);
        pattern[i] = point;
    }
    double* data_fft = new double[N];
    double* data_freq = new double[N];

    compute_psd(data, N, data_fft, data_freq);

    cpx* buffer = new cpx[sps];
    cpx* output = new cpx[N];
    double* output_fft = new double[N];
    double* output_freq = new double[N];
    double* real = new double[N];
    double* imag = new double[N];

    unsigned int* out_pattern = new unsigned int[N];

    DEBUG_INTERFACE* interf = recovery.debug();

    std::cout<<recovery.get_sample_delay();


#if DO_WINDOW
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        if(ImGui::BeginTabBar("Main Tabs")){
        if(ImGui::BeginTabItem("Item 0")){

            std::copy_n(data, N, data_c);
            for(size_t i = 0; i < N; ++i){
                channel.operate(data_c + i);
            }
            compute_psd(data_c, N, data_fft, data_freq);

            size_t k = 0;
            for(size_t i = 0; i < N; ++i){
                // recovery.operate(data_c[i], output + i);
                size_t n = recovery.operate(data_c[i], buffer);
                for(size_t j = 0; j < n; ++j){
                    output[k] = buffer[j];
                    out_pattern[k] = constel.decision(output[k]);
                    k++;
                }
            }
            compute_psd(output, N, output_fft, output_freq);

            for(size_t i = 0; i < N; ++i){
                real[i] = output[i].real();
                imag[i] = output[i].imag();
            }

            // recovery.reset();

            ImGui::SliderFloat("Frequency Offset", &offset, 0.0, 1.0);
            ImGui::SliderFloat("Noise Voltage", &v, 0.0, 1e-1, "%.5f");
            ImGui::SliderInt("Sample Delay", &symb_delay, 0, 100);
            // ImGui::SliderFloat("Loop Alpha", &alpha, 0.0, 1.0, "%.5f");
            // ImGui::SliderFloat("Loop Beta", &beta, 0.0, 1.0, "%.5f");
            // ImGui::SliderFloat("Frequency Range", &max_freq, 0.0, 1.0, "%.5f");
            channel.set_noise(v);
            channel.set_offset(offset);
            // recovery.set_pll_alpha(alpha);
            // recovery.set_pll_beta(beta);
            // recovery.set_max_freq(max_freq);
            // recovery.set_min_freq(-max_freq);

            if(ImPlot::BeginSubplots("Data", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Spectrum")){
                ImPlot::SetupAxesLimits(-0.5, 0.5, -30, 5);
                ImPlot::PlotLine("Input", data_freq, data_fft, N);
                ImPlot::PlotLine("Output", output_freq, output_fft, N);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Time Data")){
                ImPlot::SetupAxesLimits(0, points, -1.0 * static_cast<double>(constel.get_size() - 1), constel.get_size());
                ImPlot::PlotLine("Input", pattern, points);
                ImPlot::PlotLine("Output", out_pattern + symb_delay, points - symb_delay);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            }
            if(ImPlot::BeginSubplots("Debug", 1, 2, ImVec2(-1, 750))){
            if(ImPlot::BeginPlot("Constellation")){
                ImPlot::SetupAxesLimits(-2, 2, -2, 2);
                ImPlot::PlotScatter("", real + symb_delay, imag + symb_delay, N - symb_delay);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Errors")){
                ImPlot::SetupAxesLimits(0, interf->n, -2, 2);
                ImPlot::PlotLine("PLL Error", interf->PLL_ERR_HIST, interf->n);
                ImPlot::PlotLine("PLL Phase", interf->PLL_PHASE_HIST, interf->n);
                ImPlot::PlotLine("PLL Frequency", interf->PLL_FREQ_HIST, interf->n);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 1")){
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Plot 1", ImVec2(-1, 750))){
                ImPlot::EndPlot();
            }
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

