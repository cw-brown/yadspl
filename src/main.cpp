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
    size_t sps = 8;
    size_t n_filt = 32;
    size_t points = 250;
    size_t N = sps * points;

    static float v = 0.0f;
    static float offset = 0.0f;
    static int arm = 0;

    constellation_qpsk constel{};
    noise<double> sig_gen{};
    channel_model channel(offset, v);
    rectangular_modulator modulator(&constel, sps, n_filt, 0.35);
    firefighter recovery(&constel, sps, n_filt, 2.0*3.1415/100.0, 0.35);
    cpx* data = new cpx[N];
    cpx* data_c = new cpx[N];

    auto bank = recovery.get_bank();
    auto dbank = recovery.get_d_bank();

    for(size_t i = 0; i < points; ++i){
        auto point = sig_gen.random_int_range(0, constel.get_size() - 1);
        modulator.operate(point, data + i * sps);
    }

    double* data_fft = new double[N];
    double* data_freq = new double[N];

    cpx* buffer = new cpx;
    cpx* output = new cpx[points];
    double* output_fft = new double[points];
    double* output_freq = new double[points];
    double* real = new double[points];
    double* imag = new double[points];

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

            size_t j = 0;
            for(size_t i = 0; i < N; ++i){
                if(recovery.operate(data_c[i], buffer) != 0){
                    output[j] = *buffer;
                    j++;
                } 
            }
            compute_psd(output, points, output_fft, output_freq);

            for(size_t i = 0; i < points; ++i){
                real[i] = output[i].real();
                imag[i] = output[i].imag();
            }

            ImGui::SliderFloat("Frequency Offset", &offset, 0.0, 1.0);
            ImGui::SliderFloat("Noise Voltage", &v, 0.0, 1e-1, "%.5f");
            channel.set_noise(v);
            channel.set_offset(offset);

            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                ImPlot::SetupAxesLimits(-0.5, 0.5, -30, 5);
                ImPlot::PlotLine("Input", data_freq, data_fft, N);
                ImPlot::PlotLine("Output", output_freq, output_fft, points);
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Plot 1", ImVec2(-1, 750))){
                ImPlot::SetupAxesLimits(-2, 2, -2, 2);
                ImPlot::PlotScatter("", real, imag, N);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Items 1")){
            ImGui::SliderInt("Arm", &arm, 0, 31);
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                ImPlot::PlotBars("", bank[arm].get_taps(), bank[arm].get_num_taps());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Plot 1", ImVec2(-1, 750))){
                ImPlot::PlotBars("", dbank[arm].get_taps(), dbank[arm].get_num_taps());
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

