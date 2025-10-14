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
// #include "iir_filter.hpp"
// #include "polynomial.hpp"
#include "noise.hpp"
#include "yadpsl_math.hpp"
#include "constellations.hpp"
#include "symbol_rec.hpp"
#include "ring.hpp"
#include "transmitter.hpp"
#include "polyphase.hpp"
#include "graphics.hpp"
// #include "filter.hpp"
#include "modem.hpp"

#include "fft.hpp"

#define DO_WINDOW true

int main(){
    size_t sps = 4;
    size_t n_filt = 32;

    noise<double> sig_gen{};
    constellation_qpsk QPSK{};
    merm modulator(&QPSK, 15, sps, n_filt);
    auto prototype = root_nyquist(n_filt, n_filt, 1.0, 0.35, 8*sps*n_filt);
    resampler<std::complex<double>> arb(sps, n_filt, prototype);

    std::complex<double>* buffer = new std::complex<double>[sps * 256];
    std::vector<std::complex<double>> data;
    // fill a buffer completely with modulated symbols
    for(size_t i = 0; i < 256; ++i){
        // auto point = QPSK.get_point(sig_gen.randomValue(QPSK.get_bps()));
        modulator.operate(sig_gen.randomIntRange(0, QPSK.get_size() - 1), buffer + (i * sps));
    }

    double* psd = compute_psd(buffer, 256 * sps);


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
            if(ImPlot::BeginPlot("Plot 0", ImVec2(-1, 750))){
                // ImPlot::PlotLine("", y.data(), y.size());
                ImPlot::PlotLine("", psd, sps*256);
                // ImPlot::PlotStems("", h.data(), h.size());
                // ImPlot::PlotLine("", dat, 150 * sps);
                // ImPlot::PlotLine("", out, 150*sps);
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

