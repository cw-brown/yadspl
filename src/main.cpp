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
    size_t points = 256;

    constellation_qpsk QPSK{};
    noise<double> sig_gen{};
    fast_noise awgn{};
    std::complex<double>* buffer = new std::complex<double>[sps * points];
    double* fft_buffer = nullptr;
    rectangular_modulator modulator(&QPSK, sps, n_filt, 0.35);
    int k = 0;

    firefighter rec(&QPSK, sps, n_filt, 2.0*3.141/100.0, 0.35);
    rec.set_fll_size(64);
    std::complex<double>* fll_buffer = new std::complex<double>[sps * points];
    double* fll_fft = nullptr;


#if DO_WINDOW
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        static float v = 0.0f;

        auto point = sig_gen.randomIntRange(0, QPSK.get_size() - 1);
        modulator.operate(point, buffer + k);
        for(size_t i = 0; i < sps; ++i){
            buffer[k + i] += awgn.noise_voltage(v);
        }
        fft_buffer = compute_psd(buffer, sps * points);

        for(size_t i = 0; i < sps; ++i){
            rec.operate(buffer[k + i], fll_buffer + (k + i));
        }
        fll_fft = compute_psd(fll_buffer, sps * points);


        k = (k + sps) % (sps * points);

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        if(ImGui::BeginTabBar("Main Tabs")){
        if(ImGui::BeginTabItem("Item 0")){
            ImGui::SliderFloat("Noise Power", &v, 0.0f, 1.0f);
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                ImPlot::PlotLine("TX Stream FFT", fft_buffer, sps*points);
                ImPlot::PlotLine("FLL FFT", fll_fft, sps*points);
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

