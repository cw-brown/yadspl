#include <iostream>
#include <complex>
#include <array>
#include <memory>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
#include "polynomial.hpp"
#include "noise.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

int main(){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    int sps = 4;
    int span = 12;
    float beta = 0.5f;

    float xaxis[5000];

    fir_filter<float, filter_type::Root_Raised_Cosine> filt;
    for(size_t i = 0; i < 5000; ++i)
        xaxis[i] = i;


    while(!glfwWindowShouldClose(window)){
        glfwPollEvents();
        glfw_frame();

        filt.design(sps, beta, span);
        filt.response(5000);


        ImGui::Begin("Filter Impulse Response Plot");
        ImGui::SliderInt("Samples Per Symbol", &sps, 2, 16);
        ImGui::SliderInt("Filter Span", &span, 1, 32);
        ImGui::SliderFloat("Filter Roll-off", &beta, 0.001f, 1.0f);
        if(ImPlot::BeginPlot("Filter Response")){
            ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxesLimits(0, filt.getNumTaps(), -1, 1);
            ImPlot::PlotBars("Root Raised Cosine", filt.getTaps(), filt.getNumTaps());
            ImPlot::EndPlot();
        }
        if(ImPlot::BeginPlot("Frequency Response")){
            ImPlot::SetupAxes("Frequency (radians)", "Amplitude (dB)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxesLimits(0, 3.141592, -80, 10);
            ImPlot::PlotLine("Root Raised Cosine", filt.getFrequencyPoints(), filt.getMagnitudeResponse(), 5000);
            ImPlot::EndPlot();
        }
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();

    return 0;
}

