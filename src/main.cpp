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

#define DO_WINDOW true

int main(){
    size_t eval_points = 5000;

    fir_filter<float, filter_type::Lowpass_Ideal> lowpass;
    int lowpass_order = 21;
    float lowpass_fc = 0.5;

    fir_filter<float, filter_type::Highpass_Ideal> highpass;
    int highpass_order = 21;
    float highpass_fc = 0.5;

    fir_filter<float, filter_type::Bandpass_Ideal> bandpass;
    int bandpass_order = 21;
    float bandpass_fc1 = 0.25;
    float bandpass_fc2 = 0.75;

    fir_filter<float, filter_type::Bandstop_Ideal> bandstop;
    int bandstop_order = 21;
    float bandstop_fc1 = 0.25;
    float bandstop_fc2 = 0.75;

    fir_filter<float, filter_type::Lowpass_Ideal, window_type::Tukey> filt;
    int test_order = 21;
    float test_fc = 0.5;

    // filt.design(test_order, test_fc);
    // filt.synthesize();

    // for(size_t i = 0; i < filt.getNumTaps(); ++i){
    //     std::cout<<filt.getTaps()[i]<<" ";
    // }

    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);


    while(!glfwWindowShouldClose(window)){
        glfwPollEvents();
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;


        ImGui::Begin("Plottings", (bool*)true, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("Testing")){
            filt.design(test_order, test_fc);
            filt.synthesize();
            filt.response(eval_points);
            ImGui::SliderInt("Filter Order", &test_order, 2, 250);
            ImGui::SliderFloat("Cutoff Frequency", &test_fc, 0.0f, 3.14159265f);
            if(ImPlot::BeginPlot("Impulse Response")){
                ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotBars("Windowed Lowpass", filt.getTaps(), filt.getNumTaps());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Frequency Response")){
                ImPlot::SetupAxes("Frequency (Radians)", "Magnitude (dB)", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxesLimits(0.0, 3.141592, -80, 10);
                ImPlot::PlotLine("Windowed Lowpass", filt.getFrequencyPoints(), filt.getMagnitudeResponse(), eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Lowpass Filter")){
            lowpass.design(lowpass_order, lowpass_fc);
            lowpass.response(eval_points);
            ImGui::SliderInt("Filter Order", &lowpass_order, 2, 250);
            ImGui::SliderFloat("Cutoff Frequency", &lowpass_fc, 0.0f, 3.14159265f);
            if(ImPlot::BeginPlot("Impulse Response")){
                ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotBars("Lowpass", lowpass.getTaps(), lowpass.getNumTaps());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Frequency Response")){
                ImPlot::SetupAxes("Frequency (Radians)", "Magnitude (dB)", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxesLimits(0.0, 3.141592, -80, 10);
                ImPlot::PlotLine("Lowpass", lowpass.getFrequencyPoints(), lowpass.getMagnitudeResponse(), eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Highpass Filter")){
            highpass.design(highpass_order, highpass_fc);
            highpass.response(eval_points);
            ImGui::SliderInt("Filter Order", &highpass_order, 2, 250);
            ImGui::SliderFloat("Cutoff Frequency", &highpass_fc, 0.0f, 3.14159265f);
            if(ImPlot::BeginPlot("Impulse Response")){
                ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotBars("Highpass", highpass.getTaps(), highpass.getNumTaps());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Frequency Response")){
                ImPlot::SetupAxes("Frequency (Radians)", "Magnitude (dB)", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxesLimits(0.0, 3.141592, -80, 10);
                ImPlot::PlotLine("Highpass", highpass.getFrequencyPoints(), highpass.getMagnitudeResponse(), eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Bandpass Filter")){
            bandpass.design(bandpass_order, bandpass_fc1, bandpass_fc2);
            bandpass.response(eval_points);
            ImGui::SliderInt("Filter Order", &bandpass_order, 2, 250);
            ImGui::SliderFloat("Start Frequency", &bandpass_fc1, 0.0f, 3.14159265f);
            ImGui::SliderFloat("Stop Frequency", &bandpass_fc2, 0.0f, 3.14159265f);
            if(ImPlot::BeginPlot("Impulse Response")){
                ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotBars("Bandpass", bandpass.getTaps(), bandpass.getNumTaps());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Frequency Response")){
                ImPlot::SetupAxes("Frequency (Radians)", "Magnitude (dB)", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxesLimits(0.0, 3.141592, -80, 10);
                ImPlot::PlotLine("Bandpass", bandpass.getFrequencyPoints(), bandpass.getMagnitudeResponse(), eval_points);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if(ImGui::BeginTabItem("Bandstop Filter")){
            bandstop.design(bandstop_order, bandstop_fc1, bandstop_fc2);
            bandstop.response(eval_points);
            ImGui::SliderInt("Filter Order", &bandstop_order, 2, 250);
            ImGui::SliderFloat("Start Frequency", &bandstop_fc1, 0.0f, 3.14159265f);
            ImGui::SliderFloat("Stop Frequency", &bandstop_fc2, 0.0f, 3.14159265f);
            if(ImPlot::BeginPlot("Impulse Response")){
                ImPlot::SetupAxes("Sample Number", "Amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotBars("Bandstop", bandstop.getTaps(), bandstop.getNumTaps());
                ImPlot::EndPlot();
            }
            if(ImPlot::BeginPlot("Frequency Response")){
                ImPlot::SetupAxes("Frequency (Radians)", "Magnitude (dB)", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxesLimits(0.0, 3.141592, -80, 10);
                ImPlot::PlotLine("Bandstop", bandstop.getFrequencyPoints(), bandstop.getMagnitudeResponse(), eval_points);
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

