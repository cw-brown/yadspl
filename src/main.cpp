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

// #include "filter.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW true

int main(){
    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    constellation_16qam qam{};
    auto qam16_IQ = qam.get_IQ();

    constellation_qpsk qpsk{};
    auto qpsk_IQ = qpsk.get_IQ();

    constellation_bpsk bpsk{};
    auto bpsk_IQ = bpsk.get_IQ();

    while(!glfwWindowShouldClose(window)){
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;

        ImGui::Begin("Plottings", nullptr, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("Test Feature")){
            if(ImPlot::BeginPlot("IQ", ImVec2(-1, 700))){
                ImPlot::SetupAxisLimits(ImAxis_X1, -2, 2);
                ImPlot::SetupAxisLimits(ImAxis_Y1, -2, 2);
                ImPlot::PlotScatter("16QAM", qam16_IQ.first, qam16_IQ.second, qam.get_size());
                ImPlot::PlotScatter("QPSK", qpsk_IQ.first, qpsk_IQ.second, qpsk.get_size());
                ImPlot::PlotScatter("BPSK", bpsk_IQ.first, bpsk_IQ.second, bpsk.get_size());
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

