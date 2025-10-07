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

#include "fft.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW true

int main(){
    size_t sps = 4;
    noise<double> sigGen{};
    std::vector<double> prot = root_nyquist(32, 32, 1.0, 0.35, 8*32*sps);
    polyphase_upsampler samp(sps, prot, 32);

    constellation_bpsk QPSK{};

    std::vector<std::complex<double>> data(1024);  
    for(size_t i = 0; i < 1024; ++i){
        data[i] = QPSK.get_point(sigGen.randomValue(QPSK.get_bps()));
    }

    // Transmitted data, upsampled and filtered
    std::vector<std::complex<double>> c = samp.filterN(data);
    std::vector<std::complex<double>> txdat(c);
    fft_transform(c, false);


    // FLL loop recovery of transmitted data in frequency
    // frequency_recovery freq_rec(sps, 0.35, 2.0*3.14/100.0, 64);
    // std::vector<std::complex<double>> freq_rec_dat = freq_rec.operate(txdat);
    // std::vector<std::complex<double>> d = freq_rec_dat;
    // fft_transform(d, false);

    auto freq_rec_dat = txdat;

    // PLL loop recovery of FLL recovery for phase data
    phase_recovery phase_rec(sps, 2.0*3.14/100.0, 32, 1.5, QPSK, 0.35);
    std::cout<<"Kp: "<<phase_rec.get_p_gain()<<"\n";
    std::cout<<"Ki: "<<phase_rec.get_i_gain()<<"\n";
    std::vector<std::complex<double>> phase_rec_dat = phase_rec.operate(freq_rec_dat);
    std::vector<std::complex<double>> e(phase_rec_dat);
    fft_transform(e, false);
    std::cout<<"Operation Output Size: "<<phase_rec_dat.size()<<"\n";

    // Recover some symbols (maybe)
    // std::vector<std::complex<double>> recovered_symbs(phase_rec_dat.size());
    std::vector<double> real_rec(phase_rec_dat.size());
    std::vector<double> imag_rec(phase_rec_dat.size());
    for(size_t i = 0; i < phase_rec_dat.size(); ++i){
        real_rec[i] = phase_rec_dat[i].real();
        imag_rec[i] = phase_rec_dat[i].imag();
    }

    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    static std::vector<double> plot = make_psd(c);
    // static std::vector<double> plotrx = make_psd(d, false);
    static std::vector<double> plot_phase = make_psd(e);

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize 
            | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;
    
        ImGui::Begin("Plottings", nullptr, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("TX Data")){
            if(ImPlot::BeginPlot("TX Data", ImVec2(-1, 750))){
                ImPlot::PlotLine("TX", plot.data(), plot.size());
                // ImPlot::PlotLine("RX", plotrx.data(), plotrx.size());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("RX Data", ImVec2(-1, 750))){
                ImPlot::PlotLine("Phase Recovery", plot_phase.data(), plot_phase.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("RX Constellation")){
            if(ImPlot::BeginPlot("RX", ImVec2(-1, 750))){
                ImPlot::PlotScatter("BPSK", real_rec.data(), imag_rec.data(), real_rec.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();
    }
    return 0;
}

