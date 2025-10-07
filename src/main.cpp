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
    auto prot = root_nyquist(32, 32, 1.0, 0.35, 8*32*sps);
    polyphase_upsampler samp(sps, prot, 32);

    constellation_qpsk QPSK{};

    std::vector<std::complex<double>> data(1024);  
    for(size_t i = 0; i < 1024; ++i){
        data[i] = QPSK.get_point(sigGen.randomValue(QPSK.get_bps()));
    }

    auto c = samp.filterN(data);
    auto txdat(c);
    fft_transform_radix2(c, false);
    std::vector<double> plot = make_psd(c);

    std::vector<double> real(1024);
    std::vector<double> imag(1024);
    for(size_t i = 0; i < 1024; ++i){
        real[i] = txdat[i].real();
        imag[i] = txdat[i].imag();
    }

    carrier_recovery rec(sps, 0.35, 2.0*3.14/100.0, 55);

    std::vector<std::complex<double>> rxdat = rec.operate(txdat);
    std::vector<std::complex<double>> d = rxdat;
    fft_transform_radix2(d, false);
    std::vector<double> plotrx = make_psd(d, false);

    std::vector<double> realo(rxdat.size());
    std::vector<double> imago(rxdat.size());
    for(size_t i = 0; i < rxdat.size(); ++i){
        realo[i] = rxdat[i].real();
        imago[i] = rxdat[i].imag();
    }

    std::cout<<rxdat.size();

    if(DO_WINDOW){
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
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("TX Data")){
            if(ImPlot::BeginPlot("TX Data", ImVec2(-1, 750))){
                ImPlot::PlotLine("FFT", plot.data(), plot.size());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("RX Data", ImVec2(-1, 750))){
                ImPlot::PlotLine("FFT", plotrx.data(), plotrx.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("TX Data Time")){
            if(ImPlot::BeginPlot("TX", ImVec2(-1, 750))){
                ImPlot::PlotLine("Time Real", real.data(), real.size());
                ImPlot::PlotLine("Time Imag", imag.data(), imag.size());
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("RX Data Time")){
            if(ImPlot::BeginPlot("RX", ImVec2(-1, 750))){
                ImPlot::PlotLine("Time Real", realo.data(), realo.size());
                ImPlot::PlotLine("Time Imag", imago.data(), imago.size());
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

